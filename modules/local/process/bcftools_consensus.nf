// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BCFTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::bcftools=1.11' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bcftools:1.11--h7c999a4_0'
    } else {
        container 'quay.io/biocontainers/bcftools:1.11--h7c999a4_0'
    }

    input:
    tuple val(meta), path(vcf)
    tuple val(meta), path(tbi)
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path('*.fa'), emit: fasta
    path  '*.version.txt'        , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bcftools consensus -f $fasta $vcf -o ${prefix}.consensus.fa
    bedtools genomecov \\
        -bga \\
        -ibam $bam \\
        | awk '\$4 < $params.min_coverage' | bedtools merge > ${prefix}.mask.bed
    bedtools maskfasta \\
        -fi ${prefix}.consensus.fa \\
        -bed ${prefix}.mask.bed \\
        -fo ${prefix}.consensus.masked.fa
    rename_headers_masked_consensus.py -m ${prefix}.consensus.masked.fa -p ${prefix}
    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """
}