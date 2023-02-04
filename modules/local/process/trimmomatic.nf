// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::trimmomatic=0.39' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2'
    } else {
        container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.stats')        , emit: stats
    tuple val(meta), path('*.log')          , emit: log
    path '*.version.txt'                    , emit: version


    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"


    if (params.single_end && !params.adapters) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        trimmomatic SE \\
            -threads $task.cpus \\
            -summary ${prefix}.stats \\
            ${prefix}.fastq.gz \\
            ${prefix}.trim.fastq.gz \\
            LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.window_size}:${params.window_quality} MINLEN:${params.min_length} \\
            2> ${prefix}.trimmomatic.log

        echo \$(trimmomatic -version 2>&1) | sed -e "s/trimmomatic //g" > ${software}.version.txt
        """
    } else if (params.single_end && params.adapters) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        trimmomatic SE \\
            -threads $task.cpus \\
            -summary ${prefix}.stats \\
            ${prefix}.fastq.gz \\
            ${prefix}.trim.fastq.gz \\
            ILLUMINACLIP:${params.adapters}:${params.illumina_clip}\\
            LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.window_size}:${params.window_quality} MINLEN:${params.min_length} \\
            2> ${prefix}.trimmomatic.log

        echo \$(trimmomatic -version 2>&1) | sed -e "s/trimmomatic //g" > ${software}.version.txt
        """
    } else if (!params.single_end && !params.adapters) {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trimmomatic PE \\
            -threads $task.cpus \\
            -summary ${prefix}.stats \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz \\
            ${prefix}_1.paired.trim.fastq.gz \\
            ${prefix}_1.unpaired.fastq.gz \\
            ${prefix}_2.paired.trim.fastq.gz \\
            ${prefix}_2.unpaired.fastq.gz \\
            LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.window_size}:${params.window_quality} MINLEN:${params.min_length} \\
            2> ${prefix}.trimmomatic.log

        echo \$(trimmomatic -version 2>&1) | sed -e "s/trimmomatic //g" > ${software}.version.txt
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trimmomatic PE \\
            -threads $task.cpus \\
            -summary ${prefix}.stats \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz \\
            ${prefix}_1.paired.trim.fastq.gz \\
            ${prefix}_1.unpaired.fastq.gz \\
            ${prefix}_2.paired.trim.fastq.gz \\
            ${prefix}_2.unpaired.fastq.gz \\
            ILLUMINACLIP:${params.adapters}:${params.illumina_clip} \\
            LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.window_size}:${params.window_quality} MINLEN:${params.min_length} \\
            2> ${prefix}.trimmomatic.log

        echo \$(trimmomatic -version 2>&1) | sed -e "s/trimmomatic //g" > ${software}.version.txt
        """
    }
}