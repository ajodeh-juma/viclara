/*
 * Mark duplicates with Picard Markduplicates
 */

params.picard_markduplicates_options    = [:]
params.samtools_options                 = [:]

include { PICARD_MARKDUPLICATES   }       from '../../nf-core/software/picard/markduplicates/main'           addParams( options: params.picard_markduplicates_options    )
include { BAM_SORT_SAMTOOLS       }       from '../../nf-core/subworkflow/bam_sort_samtools'                 addParams( options: params.samtools_options )

workflow MARKDUPLICATES_PICARD {
    take:
    bam       // channel: [ val(meta), [ bam ] ]
    
    main:
    /*
     * Mark duplicates with PICARD MARKDUPLICATES
     */
    PICARD_MARKDUPLICATES ( bam )

    /*
     * Run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( PICARD_MARKDUPLICATES.out.bam )

    emit:
    bam              = BAM_SORT_SAMTOOLS.out.bam           // channel: [ val(meta), bam   ]
    metrics          = PICARD_MARKDUPLICATES.out.metrics    // channel: [ val(meta), metrics ]
    picard_version   = PICARD_MARKDUPLICATES.out.version       // path: *.version.txt

    bai              = BAM_SORT_SAMTOOLS.out.bai       // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats     // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat  // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats  // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version   // path: *.version.txt
}