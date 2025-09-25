#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/viclara
========================================================================================
 nf-core/viclara Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/viclara
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl = 2

//log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run main.nf --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}


////////////////////////////////////////////////////
/* --        GENOME PARAMETER VALUES           -- */
////////////////////////////////////////////////////

// params.fasta         = Checks.get_genome_attribute(params, 'fasta')
// params.bed           = Checks.get_genome_attribute(params, 'bed')
// params.gtf           = Checks.get_genome_attribute(params, 'gtf')
// params.gff           = Checks.get_genome_attribute(params, 'gff')
// params.bwa_index     = Checks.get_genome_attribute(params, 'bwa')
// params.bowtie2_index = Checks.get_genome_attribute(params, 'bowtie2')


////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

// Check genome key exists if provided
// Checks.genome_exists(params, log)


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

include { viclara } from './viclara' addParams( summary_params: summary_params )

workflow {
    viclara ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+

// if (params.validate_params) {
//     NfcoreSchema.validateParameters(params, json_schema, log)
// }
