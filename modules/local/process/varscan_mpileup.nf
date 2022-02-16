// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process VARSCAN_MPILEUP {
    label 'process_medium'
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::varscan=2.4.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/varscan:2.4.4--h9402c20_2"
    } else {
        container "quay.io/biocontainers/varscan:2.4.4--h9402c20_2"
    }

    input:
    tuple val(meta), path(mpileup)

    output:
    tuple val(meta), path("*.vcf.gz")           , emit: vcf
    tuple val(meta), path("*.tbi")              , emit: tbi
    tuple val(meta), path("*stats.txt")         , emit: stats
    tuple val(meta), path("*.varscan.log")      , emit: log

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def lofreq_prefix = "${prefix}.AF${params.min_allele_freq}"
    def hifreq_prefix = "${prefix}.AF${params.max_allele_freq}"
    def strand = params.varscan_strand_filter ? "--strand-filter 1" : "--strand-filter 0"

    """
    echo "${meta.id}" > sample_name.list
    varscan mpileup2cns \\
        $mpileup \\
        $options.args \\
        --output-vcf 1 \\
        --variants \\
        --vcf-sample-list sample_name.list \\
        $strand \\
        2> ${prefix}.varscan.log \\
        | bgzip -c > ${lofreq_prefix}.vcf.gz
    tabix -p vcf -f ${lofreq_prefix}.vcf.gz
    bcftools stats ${lofreq_prefix}.vcf.gz > ${lofreq_prefix}.bcftools_stats.txt
    sed -i.bak '/LC_ALL/d' ${prefix}.varscan.log

    bcftools filter \\
        $options.args2 \\
        --output-type z \\
        --output ${hifreq_prefix}.vcf.gz \\
        ${lofreq_prefix}.vcf.gz
    tabix -p vcf -f ${hifreq_prefix}.vcf.gz
    bcftools stats ${hifreq_prefix}.vcf.gz > ${hifreq_prefix}.bcftools_stats.txt
    varscan 2>&1 | sed 's/^.*VarScan //; s/ .*\$//' | sed -n 1p > ${software}.version.txt
    """
}