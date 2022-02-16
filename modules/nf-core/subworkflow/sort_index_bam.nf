/*
 * Sort, index BAM file and run samtools stats, flagstat and idxstats
 */

params.options = [:]

include { SAMTOOLS_SORT      } from '../software/samtools/sort/main'  addParams( options: params.options )
include { SAMTOOLS_INDEX     } from '../software/samtools/index/main' addParams( options: params.options )

workflow SORT_INDEX_BAM {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    
    main:
    SAMTOOLS_SORT      ( ch_bam )
    SAMTOOLS_INDEX     ( SAMTOOLS_SORT.out.bam )
    
    emit:
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    version  = SAMTOOLS_SORT.out.version       // path: *.version.txt
}