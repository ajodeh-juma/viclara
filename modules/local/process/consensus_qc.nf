// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options) 

process CONSENSUS_QC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::python=3.6.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.6.1"
    } else {
        container "quay.io/biocontainers/python:3.6.1"
    }

    input:
    tuple val(meta), path(bam), path(fasta)
    path reference

    output:
    tuple val(meta), path("*.qc.csv")   , emit: csv


    script:
    def prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    generate_consensus_qc.py \\
        --prefix $prefix \\
        --bam $bam \\
        --reference $reference \\
        --consensus-fasta $fasta \\
        --outfile ${prefix}.qc.csv
    """
}