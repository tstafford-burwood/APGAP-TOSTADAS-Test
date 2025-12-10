#!/usr/bin/env nextflow 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                    RUNNING SUBMISSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREP_SUBMISSION       } from '../../modules/local/prep_submission/main'
include { SUBMIT_SUBMISSION     } from '../../modules/local/submit_submission/main'

workflow SUBMISSION {
    take:
        submission_ch         // (meta: [batch_id], samples: [ [meta, fq1, fq2, nnp], ... ]), enabledDatabases (list), batch_tsv (file)
        submission_config

    main:
        submission_config_file = file(submission_config)

        // Split submission_ch into separate channels for PREP_SUBMISSION inputs
        prep_tuple_ch = submission_ch.map { meta, samples, enabledDatabases, batch_tsv -> 
            tuple(meta, samples, enabledDatabases)
        }
        prep_batch_tsv_ch = submission_ch.map { meta, samples, enabledDatabases, batch_tsv -> 
            batch_tsv
        }

        PREP_SUBMISSION(prep_tuple_ch, prep_batch_tsv_ch, submission_config_file)

        PREP_SUBMISSION.out.submission_files
            .set { submission_batch_folder }

        SUBMIT_SUBMISSION(submission_batch_folder, submission_config_file)

        // Print message if dry_run is enabled
        if (params.dry_run) {
            log.info "Workflow ends here because dry_run is set to true. Please see submission log files for details."
        }

    emit:
        submission_batch_folder = SUBMIT_SUBMISSION.out.submission_batch_folder
}
