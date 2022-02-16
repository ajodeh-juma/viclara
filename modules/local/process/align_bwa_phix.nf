/*
 * Alignment with BWA
 */

params.bwa_mem_options    = [:]
params.sort_index_options = [:]

include { BWA_MEM           }       from '../../nf-core/software/bwa/mem/main'           addParams( options: params.bwa_mem_options    )
include { SORT_INDEX_BAM }       from '../../nf-core/subworkflow/sort_index_bam'           addParams( options: params.sort_index_options )

workflow ALIGN_BWA_PHIX {
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
    SORT_INDEX_BAM ( BWA_MEM.out.bam )

    emit:
    bam              = SORT_INDEX_BAM.out.bam
    bai              = SORT_INDEX_BAM.out.bai         //    channel: [ val(meta), [ bai ] ]
    bwa_version      = BWA_MEM.out.version            //    path: *.version.txt
    samtools_version = SORT_INDEX_BAM.out.version     //    path: *.version.txt
}