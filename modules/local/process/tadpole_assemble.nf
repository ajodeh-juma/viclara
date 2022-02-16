// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TADPOLE_ASSEMBLE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::bbmap=38.87' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bbmap:38.87--hf29c6f4_0"
    } else {
        container "quay.io/biocontainers/bbmap:38.87--hf29c6f4_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.contigs.fasta")   , emit: contigs
    //path  '*.version.txt'                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads    = params.single_end ? "-in=$reads" : "-in=${reads[0]} -in2=${reads[1]}"
    def mem = 4
    if (!task.memory) {
        log.info "[tadpole.sh] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this."
    } else {
        mem = task.memory.giga
    }

    """
    tadpole.sh \\
        threads=${task.cpus} \\
        showstats=t \\
        k=31 \\
        prealloc=t \\
        overwrite=f \\
        -Xmx${mem}g \\
        $input_reads \\
        out=${prefix}.contigs.fasta \\
        mode=contig \\
        mincountseed=3 \\
        mincountextend=2 \\ 
        minextension=2 \\
        mincontig=auto \\
        mincoverage=3 \\
        ecc=t \\
        reassemble=t
    """
}
