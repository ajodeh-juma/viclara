// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process COUNT_READS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options) }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.tsv')         , emit: tsv

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (params.single_end) {
        """
        count_reads.py --fwd ${reads[0]} --output ${prefix}.tsv
        """
    } else {
        """
        count_reads.py --fwd ${reads[0]} --rev ${reads[1]} --output ${prefix}.tsv
        """
    }
}