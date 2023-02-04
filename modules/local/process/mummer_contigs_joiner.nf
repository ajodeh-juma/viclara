// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MUMMER_CONTIGS_JOINER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::mummer=3.23' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mummer:3.23--h6de7cb9_11"
    } else {
        container "quay.io/biocontainers/mummer:3.23--h6de7cb9_11"
    }
    
    input:
    tuple val(meta), path(contigs)
    path fasta

    output:
    tuple val(meta), path("*.delta")                           , emit: delta
    tuple val(meta), path("*.pseudo.fasta")                    , emit: pseudo
    tuple val(meta), path("*.tiling.txt")                      , emit: tiling
    path '*.version.txt'                    , emit: version
    tuple val(meta), path("*.gap.filled.assembly.gaps.bed")    , emit: bed
    tuple val(meta), path("*.gap.filled.assembly.fasta")       , emit: fasta


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    nucmer \\
        -prefix ${prefix} \\
        $fasta \\
        $contigs

    show-tiling \\
        -p ${prefix}.pseudo.fasta \\
        ${prefix}.delta > ${prefix}.tiling.txt
    
    contigs_joiner.py \\
        --reference ${fasta} \\
        --sample ${prefix} \\
        --tiling ${prefix}.tiling.txt \\
        --contigs ${contigs}

    echo \$(nucmer --version 2>&1) | sed 's/^.*NUCmer (NUCleotide MUMmer) version v//; s/ .*\$//' > ${software}.version.txt
    """

}