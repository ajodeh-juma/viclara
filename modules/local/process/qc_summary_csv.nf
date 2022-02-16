// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QC_SUMMARY_CSV {
    tag "${prefix}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options) }

    input:
    path csv
    path readcounts
    val  prefix

    output:
    path  '*.csv' , emit: csv

    script:
    """
    summarize_qc_csv.py --input $csv --read-counts $readcounts --prefix $prefix
    """
}