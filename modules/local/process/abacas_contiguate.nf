// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ABACAS_CONTIGUATE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::abacas=1.3.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/abacas:1.3.1--pl526_0"
    } else {
        container "quay.io/biocontainers/abacas:1.3.1--pl526_0"
    }

    input:
    tuple val(meta), path(contigs)
    path  fasta

    output:
    tuple val(meta), path('*.nucmer.delta*')      , emit: nucmer_delta
    tuple val(meta), path('*.filtered.delta*')    , emit: filtered_delta
    tuple val(meta), path('*.tiling')             , emit: tiling
    tuple val(meta), path('*.out')                , emit: out
    tuple val(meta), path('*.bin')                , emit: bin
    tuple val(meta), path('*.fasta')              , emit: fasta
    tuple val(meta), path('*.gaps')               , emit: gaps
    tuple val(meta), path('*.gaps.tab*')          , emit: gaps_tab
    path '*.version.txt'                          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    abacas.pl \\
        -r $fasta \\
        -q $contigs \\
        $options.args \\
        -o ${prefix}.abacas
    mv nucmer.delta ${prefix}.abacas.nucmer.delta
    mv nucmer.filtered.delta ${prefix}.abacas.nucmer.filtered.delta
    mv nucmer.tiling ${prefix}.abacas.nucmer.tiling
    mv unused_contigs.out ${prefix}.abacas.unused.contigs.out
    rename_headers_masked_consensus.py -m ${prefix}.abacas.fasta -p ${prefix}
    echo \$(abacas.pl -v 2>&1) | sed 's/^.*ABACAS.//; s/ .*\$//' > ${software}.version.txt
    """
}