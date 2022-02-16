// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MERGE_COUNTS {
    tag "${prefix}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options) }

    input:
    path tsv
    val  prefix

    output:
    path  '*.csv' , emit: csv

    script:
    
    """
    merge_read_counts.py --input $tsv --prefix $prefix $options.args 
    """
}