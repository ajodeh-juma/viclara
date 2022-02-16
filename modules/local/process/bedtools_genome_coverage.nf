// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BEDTOOLS_GENOME_COVERAGE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::bedtools=2.30.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.genomecov.coverage.bed')        , emit: coverage_bed
    path  '*.version.txt'                                    , emit: version


    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    echo -e \"reference_name\tposition\tdepth\" > ${prefix}.genomecov.coverage.bed
    bedtools genomecov -d -ibam $bam >> ${prefix}.genomecov.coverage.bed
    
    bedtools --version | sed  -e 's/bedtools v//g' > ${software}.version.txt
    """
}