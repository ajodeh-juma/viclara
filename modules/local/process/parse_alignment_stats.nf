// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PARSE_FLAGSTAT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options) }

    input:
    tuple val(meta), path(flagstat)

    output:
    tuple val(meta), path('*.tsv')         , emit: tsv

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    extract_alignment_stats.py --flagstat $flagstat --output ${prefix}.tsv
    """
}