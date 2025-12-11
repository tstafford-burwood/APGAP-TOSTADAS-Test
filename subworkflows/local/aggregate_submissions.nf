include { FETCH_REPORTS                              } from '../../modules/local/fetch_reports/main'
include { AGGREGATE_REPORTS                          } from '../../modules/local/aggregate_reports/main'
include { JOIN_ACCESSIONS_WITH_METADATA              } from '../../modules/local/join_accessions_with_metadata/main'

workflow AGGREGATE_SUBMISSIONS {
    take:
      submission_dirs // works for one or more batch_dir(s)
      submission_config
      validated_metadata_tsv
      wait_signal

    main:
      FETCH_REPORTS(submission_dirs, file(submission_config))

      // Collect the individual batch submission_reports
      FETCH_REPORTS.out.submission_report
                .collect()
                .set { all_report_csvs }

      // Log if no reports are available (expected in dry_run mode)
      all_report_csvs
        .ifEmpty {
          if (params.dry_run) {
            log.info "Skipping AGGREGATE_REPORTS and JOIN_ACCESSIONS_WITH_METADATA because dry_run is enabled (no reports were fetched to aggregate)."
          } else {
            log.warn "WARNING: No CSV reports found from FETCH_REPORTS. AGGREGATE_REPORTS and JOIN_ACCESSIONS_WITH_METADATA will be skipped."
          }
        }

      // Concatenate batch csvs for all samples
      // Note: If all_report_csvs is empty (e.g., in dry_run mode), this process will be automatically skipped by Nextflow
      AGGREGATE_REPORTS(all_report_csvs)

      // Concatenate the batch TSVs, then add the (optional) accession IDs to them
      // Note: If AGGREGATE_REPORTS was skipped, this will also be skipped
      JOIN_ACCESSIONS_WITH_METADATA(validated_metadata_tsv, AGGREGATE_REPORTS.out.submission_report)

    emit:
      all_report_csvs = all_report_csvs
      accession_augmented_xlsx = params.dry_run ? Channel.empty() : JOIN_ACCESSIONS_WITH_METADATA.out.updated_excel
}