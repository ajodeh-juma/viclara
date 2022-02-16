/*
 * Alignment with BOWTIE2
 */

params.align_options    = [:]
params.samtools_options = [:]

include { BOWTIE2_ALIGN           }       from '../../nf-core/software/bowtie2/align/main'           addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS       }       from '../../nf-core/subworkflow/bam_sort_samtools'         addParams( options: params.samtools_options )

workflow ALIGN_BOWTIE2 {
    take:
    reads       // channel: [ val(meta), [ reads ] ]
    index       // channel: /path/to/bwa/index/
    
    main:
    /*
     * Map reads with BOWTIE2
     */
    BOWTIE2_ALIGN ( reads, index )

    /*
     * Run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( BOWTIE2_ALIGN.out.bam )

    emit:
    bam              = BAM_SORT_SAMTOOLS.out.bam       // channel: [ val(meta), bam   ]
    bowtie2_version  = BOWTIE2_ALIGN.out.version       // path: *.version.txt
    bai              = BAM_SORT_SAMTOOLS.out.bai       // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats     // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat  // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats  // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version   // path: *.version.txt
}