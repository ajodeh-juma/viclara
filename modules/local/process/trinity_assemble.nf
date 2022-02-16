// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRINITY_ASSEMBLE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::trinity=date.2011_11_26' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/trinity:date.2011_11_26"
    } else {
        container "quay.io/biocontainers/trinity:date.2011_11_26"
    }

    input:
    tuple val(meta), path(reads)

    output:
    // tuple val(meta), path("${prefix}_trinity/Trinity.fasta")   , emit: contig
    tuple val(meta), path("*.fasta")   , emit: contig
    path  '*.version.txt'                                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def avail_mem = 8
    if (!task.memory) {
        log.info '[Trinity] Available memory not known - defaulting to 8G. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    
    
    if (params.single_end) {
        """
        Trinity \\
            --seqType fq \\
            --single $reads \\
            --CPU ${task.cpus} \\
            --max_memory ${avail_mem}G \\
            $options.args \\
            --output ${prefix}_trinity
        mv ${prefix}_trinity/Trinity.fasta ${prefix}_trinity/${prefix}.fasta
        """
    } else {
        """
        Trinity \\
            --seqType fq \\
            --left ${reads[0]} \\
            --right ${reads[1]} \\
            --CPU ${task.cpus} \\
            --max_memory ${avail_mem}G \\
            $options.args \\
            --output ${prefix}_trinity 
        
        mv ${prefix}_trinity/Trinity.fasta ${prefix}_trinity/${prefix}.fasta
        """
    }
}
