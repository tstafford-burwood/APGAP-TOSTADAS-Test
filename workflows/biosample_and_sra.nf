#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
						 GET NECESSARY MODULES OR SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// get the utility processes / subworkflows
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { METADATA_VALIDATION                               } from "../modules/local/metadata_validation/main"
include { CHECK_VALIDATION_ERRORS							} from "../modules/local/check_validation_errors/main"
include { WRITE_VALIDATED_FULL_TSV                          } from "../modules/local/write_validated_full_tsv/main"
include { SUBMISSION		                                } from "../subworkflows/local/submission"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
									MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Global method to trim white space from file paths
def trimFile(path) {
    return path?.trim() ? file(path.trim()) : null
}

workflow BIOSAMPLE_AND_SRA {
	main:
	// Validate input parameters
	validateParameters()

	// Print summary of supplied parameters
	log.info paramsSummaryLog(workflow)

	// Run metadata validation process
	METADATA_VALIDATION ( 
		file(params.meta_path),
		file(params.submission_config),
		file(params.custom_fields_file),
		file(params.biosample_fields_key)
	)

	// Enforce error checking before anything else continues
    CHECK_VALIDATION_ERRORS(METADATA_VALIDATION.out.errors)

    // Get status from the check
	CHECK_VALIDATION_ERRORS.out.status.subscribe { status ->
		if (status == "ERROR") {
			log.info "Validation failed. Please check ${params.outdir}/${params.metadata_basename}/${params.validation_outdir}/error.txt"
			workflow.abort()
		}
	}

	// Create metadata_batch_ch from validation output
	metadata_batch_ch = METADATA_VALIDATION.out.tsv_files
		.flatten()
		.map { batch_tsv ->
			def meta = [
				batch_id: batch_tsv.getBaseName(),
				batch_tsv: batch_tsv
			]
			log.info "Created metadata_batch_ch entry: batch_id=${meta.batch_id}, file=${batch_tsv}"
			[meta, batch_tsv]
		}
		
	// Log metadata_batch_ch contents
	metadata_batch_ch
		.ifEmpty {
			log.error "ERROR: metadata_batch_ch is empty. No TSV files were created by METADATA_VALIDATION."
		}
		
	// Create a separate channel for the join (recreated from source to avoid consumption issues)
	metadata_batch_ch_for_join = METADATA_VALIDATION.out.tsv_files
		.flatten()
		.map { batch_tsv ->
			def meta = [
				batch_id: batch_tsv.getBaseName()
			]
			log.info "Created metadata_batch_ch_for_join entry: batch_id=${meta.batch_id}, file=${batch_tsv}"
			[meta, batch_tsv]
		}
		
	// Log metadata_batch_ch_for_join contents
	metadata_batch_ch_for_join
		.ifEmpty {
			log.error "ERROR: metadata_batch_ch_for_join is empty. Cannot perform join operation."
		}

	// Aggregate the tsvs for concatenation
	METADATA_VALIDATION.out.tsv_files
		.collect()
		.set { validated_tsvs_list }

	WRITE_VALIDATED_FULL_TSV ( validated_tsvs_list )
		
	if (params.submission) { 
		// Generate the (per-sample) fasta and fastq paths
		sample_ch = metadata_batch_ch.flatMap { meta, batch_tsv_file -> 
			def rows = batch_tsv_file.splitCsv(header: true, sep: '\t')
			if (rows.isEmpty()) {
				log.error "ERROR: Batch ${meta.batch_id} TSV file is empty or has no data rows. File: ${batch_tsv_file}"
				return []
			}
			log.info "Processing batch ${meta.batch_id}: Found ${rows.size()} rows in TSV file"
			return rows
				.findAll { row -> 
					def sample_id = row.sample_name?.trim()
					if (!sample_id) {
						log.warn "WARNING: Skipping row with empty sample_name in batch ${meta.batch_id}"
						return false
					}
					return true
				}
				.collect { row -> 
					def sample_id = row.sample_name?.trim()
					def fq1   = row.containsKey('int_illumina_sra_file_path_1') ? trimFile(row.int_illumina_sra_file_path_1) : null
					def fq2   = row.containsKey('int_illumina_sra_file_path_2') ? trimFile(row.int_illumina_sra_file_path_2) : null
					def nnp   = row.containsKey('int_nanopore_sra_file_path_1') ? trimFile(row.int_nanopore_sra_file_path_1) : null

					def sample_meta = [
						batch_id  : meta.batch_id,
						sample_id: sample_id
					]
					return [sample_meta, fq1, fq2, nnp]
				}
		}
		
		// Log if sample_ch is empty and provide diagnostics
		sample_ch
			.ifEmpty {
				log.error "ERROR: sample_ch is empty. No samples found in metadata TSV files."
				log.error "This means either:"
				log.error "  1. The TSV files created by METADATA_VALIDATION are empty"
				log.error "  2. The TSV files don't have a 'sample_name' column"
				log.error "  3. All 'sample_name' values in the TSV files are empty/null"
				log.error "Check the batched_tsvs/*.tsv files in the validation output directory."
			}
			// Log samples found without consuming the channel
			.map { meta, fq1, fq2, nnp ->
				log.info "Sample found: ${meta.sample_id} in batch ${meta.batch_id}"
				tuple(meta, fq1, fq2, nnp)
			}

		// Check for valid submission inputs and make batch channel
		submission_batch_ch = sample_ch
			.map { meta, fq1, fq2, nnp ->
				[meta.batch_id, [meta: meta, fq1: fq1, fq2: fq2, nanopore: nnp]]
			}
			.groupTuple()
			.map { batch_id, sample_maps ->

				def enabledDatabases = [] as Set
				def sraWarnings = [] as List

				sample_maps.each { sample ->
					def sid = sample.meta.sample_id
					// Check if file exists - handle both local and GS paths
					def fq1Exists = sample.fq1 && (sample.fq1.toString().startsWith('gs://') || file(sample.fq1).exists())
					def fq2Exists = sample.fq2 && (sample.fq2.toString().startsWith('gs://') || file(sample.fq2).exists())
					def nnpExists = sample.nanopore && (sample.nanopore.toString().startsWith('gs://') || file(sample.nanopore).exists())

					def hasIllumina = fq1Exists && fq2Exists
					def hasNanopore = nnpExists

					if (params.sra && (hasIllumina || hasNanopore)) {
						enabledDatabases << "sra"
					}
					if (params.sra && !(hasIllumina || hasNanopore)) {
						sraWarnings << sid
					}
					if (params.biosample) {
						enabledDatabases << "biosample"
					}
				}

				if (sraWarnings) {
					log.warn "SRA submission will be skipped for batch ${batch_id} due to missing data for samples: ${sraWarnings.join(', ')}"
				}

				def meta = [
					batch_id : batch_id
				]

				return tuple(meta, sample_maps, enabledDatabases as List)
			}
			// Log before join and prepare for join (batch_id must be first element for join to work)
			.map { meta, sample_maps, enabledDatabases ->
				log.info "Before join - Batch ${meta.batch_id}: ${sample_maps.size()} samples, enabledDatabases: ${enabledDatabases}"
				// Put batch_id as first element for join to match on
				[meta.batch_id, meta, sample_maps, enabledDatabases]
			}
			// Join with metadata_batch_ch_for_join to get the batch_tsv file
			.join(metadata_batch_ch_for_join.map { meta, batch_tsv -> 
				log.info "Join key: batch_id=${meta.batch_id}, batch_tsv=${batch_tsv}"
				[meta.batch_id, batch_tsv] 
			})
			.map { batch_id, meta, sample_maps, enabledDatabases, batch_tsv ->
				log.info "After join - Batch ${batch_id}: ${sample_maps.size()} samples, enabledDatabases: ${enabledDatabases}, batch_tsv=${batch_tsv}"
				// Log warning if enabledDatabases is empty
				if (enabledDatabases.isEmpty()) {
					log.warn "WARNING: No databases enabled for batch ${batch_id}. Submission will be skipped. Check that params.biosample or params.sra are set, and that sample files exist."
				}
				// Reconstruct meta with batch_id
				def final_meta = [
					batch_id: batch_id
				]
				tuple(final_meta, sample_maps, enabledDatabases, batch_tsv)
			}

		// Log if submission_batch_ch is empty
		submission_batch_ch
			.ifEmpty { 
				log.error "ERROR: submission_batch_ch is empty. No batches to submit. This will cause PREP_SUBMISSION and SUBMIT_SUBMISSION to be skipped."
				log.error "Possible causes:"
				log.error "  1. sample_ch is empty (no samples in TSV files)"
				log.error "  2. Join failed (batch_ids don't match between channels)"
				log.error "  3. enabledDatabases is empty for all batches"
			}

		SUBMISSION(
			submission_batch_ch, // meta: [batch_id], samples: [ [meta, fq1, fq2, nnp], ... ]), enabledDatabases (list), batch_tsv (file)
			params.submission_config
		)
	}

	emit:
	validated_concatenated_tsv = WRITE_VALIDATED_FULL_TSV.out.validated_concatenated_tsv // contains data for all batches of samples
	submission_batch_folder = params.submission ? SUBMISSION.out.submission_batch_folder : null
}
