// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KRAKEN2_REPORTS {
    tag "KRAKEN2_MPA_REPORTS_SUMMARY"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

    conda (params.enable_conda ? 'bioconda::kraken2=2.1.1 conda-forge::pigz=2.6' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:941789bd7fe00db16531c26de8bf3c5c985242a5-0'
    }

    input:
    path reports
    

    output:
    path ("*.tables.txt")    , emit: tables
    path ("*.species.txt")   , emit: species

    script:
    def software     = getSoftwareName(task.process)
    def prefix = "merged"
    
    """
    merge_mpa_tables.py --input $reports --prefix ${prefix}
    """
}