// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PLOT_READ_COUNTS {
    tag "${prefix}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options) }
    
    input:
    path            counts
    val             prefix
    

    output:
    path "*.pdf"

    script:
    """
    plot_read_counts.R --counts $counts --prefix $prefix --outdir ./
    """

    // if (!params.run_id) {
    //     def prefix = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
    //     """
    //     plot_read_counts.R --counts $counts --prefix $prefix --outdir ./
    //     """
    // } else {
    //     def prefix   = "${params.run_id}"
    //     """
    //     plot_read_counts.R --counts $counts --prefix $prefix --outdir ./
    //     """
    // }
}