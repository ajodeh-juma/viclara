// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process SNPEFF_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::snpeff=5.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1'
    } else {
        container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'
    }

    input:
    tuple val(meta), path(vcf)
    path fasta
    path gff

    output:
    tuple val(meta), path("*.snpEff.csv")   , emit: csv
    tuple val(meta), path("*.vcf.gz*")      , emit: vcf
    tuple val(meta), path("*.txt")          , emit: txt
    tuple val(meta), path("*.html")         , emit: html

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // def lastPath   = params.fasta.lastIndexOf(File.separator)
    // def lastExt    = params.fasta.lastIndexOf(".")
    // def index_base = params.fasta.substring(lastPath+1,lastExt)

    """
    annotate_variants.py \\
        --reference $fasta \\
        --prefix ${prefix} \\
        --gff $gff \\
        --vcf-file $vcf \\
        --csv ${prefix}.snpEff.csv \\
        --vcf-out ${prefix}.snpEff.vcf.gz

    echo \$(snpEff -version 2>&1) | sed 's/^.*SnpEff //; s/ .*\$//' > ${software}.version.txt
    """
}