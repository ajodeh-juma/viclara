// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process COVERAGE_PLOTS {
    tag "${prefix}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options) }
    
    input:
    path            consensus_fasta
    path            genomecov_bed
    file            metadata
    file            consensus_summary
    val             prefix

    output:
    path "*.pdf"

    script:
    def meta   = params.metadata ? "--metadata ${metadata}" :  " "

    // if (params.metadata) {
    //     """
    //     plot_nucleotide_densities.r --consensus_files $consensus_fasta --prefix $prefix --outdir ./
    //     plot_genome_coverage.r --coverage_files $genomecov_bed $meta --prefix $prefix --outdir ./
    //     plot_coverage_vs_ct.r --summary $consensus_summary $meta --prefix $prefix --outdir ./
    //     """
    // } else {
    // """
    // plot_nucleotide_densities.r --consensus_files $consensus_fasta --prefix $prefix --outdir ./
    // plot_genome_coverage.r --coverage_files $genomecov_bed --prefix $prefix --outdir ./
    // plot_coverage_vs_ct.r --summary $consensus_summary --prefix $prefix --outdir ./
    // """

    """
    plot_nucleotide_densities.r --consensus_files $consensus_fasta --prefix $prefix --outdir ./
    plot_genome_coverage.r --coverage_files $genomecov_bed $meta --prefix $prefix --outdir ./
    plot_coverage_vs_ct.r --summary $consensus_summary $meta --prefix $prefix --outdir ./
    """
    //}
    
}