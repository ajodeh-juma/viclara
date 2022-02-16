// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SAMTOOLS_UNMAPPED_FASTQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz")  ,  emit: reads
    path  "*.version.txt"                ,  emit: version

    
    // if (params.single_end) {
    //     """
    //     samtools view -b --threads ${task.cpus} -f 4 $bam > ${prefix}.unmapped.bam
    //     samtools sort -n --threads $task.cpus -o ${prefix}.unmapped.sorted.bam -T $prefix ${prefix}.unmapped.bam
    //     bedtools bamtofastq -i ${prefix}.unmapped.sorted.bam -fq ${prefix}_R1.fastq
    //     echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    //     """
    // } else {
    //     """
    //     samtools view -b --threads ${task.cpus} -f 4 $bam > ${prefix}.unmapped.bam
    //     samtools sort -n --threads $task.cpus -o ${prefix}.unmapped.sorted.bam -T $prefix ${prefix}.unmapped.bam
    //     bedtools bamtofastq -i ${prefix}.unmapped.sorted.bam -fq ${prefix}_R1.fastq  -fq2 ${prefix}_R2.fastq
    //     echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    //     """
    // }    

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (params.single_end) {
        """
        samtools fastq $options.args -f 4--threads $task.cpus $bam  > ${prefix}.fastq.gz

        echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
        """

    } else {
        """
        samtools fastq $options.args -f 4 --threads $task.cpus -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -0 /dev/null -s /dev/null $bam

        echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
        """
    }
}