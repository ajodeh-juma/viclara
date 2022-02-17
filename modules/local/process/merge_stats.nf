// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MERGE_STATS {
    tag "${prefix}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options) }

    input:
    path raw
    path trimmed
    path mapped
    val  prefix

    output:
    path  '*.csv' , emit: csv

    script:
    
    

    if (!params.run_id) {
        """
        merge_stats.py --raw $raw  --trimmed $trimmed --mapped $mapped --prefix $prefix
        """
    } else {
        """
        merge_stats.py --raw $raw  --trimmed $trimmed --mapped $mapped --prefix $prefix
        """
    }
}