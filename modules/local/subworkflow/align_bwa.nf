/*
 * Alignment with BWA
 */

params.align_options    = [:]
params.samtools_options = [:]

include { BWA_MEM           }       from '../../nf-core/software/bwa/mem/main'           addParams( options: params.align_options    )
include { BAM_SORT_SAMTOOLS }       from '../../nf-core/subworkflow/bam_sort_samtools'   addParams( options: params.samtools_options )

workflow ALIGN_BWA {
    take:
    reads       // channel: [ val(meta), [ reads ] ]
    index       // channel: /path/to/bwa/index/
    
    main:
    /*
     * Map reads with BWA
     */
    BWA_MEM ( reads, index )

    /*
     * Run samtools stats, flagstat and idxstats
     */
    BAM_SORT_SAMTOOLS ( BWA_MEM.out.bam )

    emit:
    bam              = BAM_SORT_SAMTOOLS.out.bam
    bwa_version      = BWA_MEM.out.version            //    path: *.version.txt
    bai              = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats            = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    samtools_version = BAM_SORT_SAMTOOLS.out.version  //    path: *.version.txt
}