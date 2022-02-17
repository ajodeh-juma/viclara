// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process READ_COUNT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options) }
    
    input:
    tuple val(meta), path(reads)    

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val(meta), path("*.csv"),   emit: csv

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (params.single_end) {
        """
        reformat.sh in=$reads -Xmx4G 2> ${prefix}.stats
        parse_read_counts.py --stats ${prefix}.stats --output ${prefix}.csv
        """
    } else {
        """
        reformat.sh in=${reads[0]} in1=${reads[1]} -Xmx4G 2> ${prefix}.stats
        parse_read_counts.py --stats ${prefix}.stats --output ${prefix}.csv
        """
    }    
}