/*
 * Uncompress and prepare reference genome files
*/

params.host_genome_options         = [:]
params.host_index_options          = [:]
params.host_bwa_index_options      = [:]
params.host_bowtie2_index_options  = [:]


include { GUNZIP as GUNZIP_FASTA           }  from '../process/gunzip'                                     addParams( options: params.host_genome_options               )
include { UNTAR as UNTAR_BWA_INDEX         }  from '../process/untar'                                      addParams( options: params.host_bwa_index_options            )
include { UNTAR as UNTAR_BOWTIE2_INDEX     }  from '../process/untar'                                      addParams( options: params.host_bowtie2_index_options        )
include { BWA_INDEX                        }  from '../../nf-core/software/bwa/index/main'                 addParams( options: params.host_bwa_index_options            )
include { BOWTIE2_BUILD                    }  from '../../nf-core/software/bowtie2/build/main'             addParams( options: params.host_bowtie2_index_options        )

workflow PREPARE_HOST_GENOME {
    take:
    prepare_tool_indices // list: tools to prepare indices for

    main:

    /* Uncompress genome fasta file if required */
    ch_host_fasta = Channel.empty()
    if (params.host_fasta) {
        if (params.host_fasta.endsWith('.gz')) {
            ch_host_fasta = GUNZIP_FASTA ( params.host_fasta ).gunzip
        } else {
            ch_host_fasta = file(params.host_fasta)
        }
    }

    /* Uncompress BOWTIE2 index or generate from scratch if required */ 
    ch_host_bowtie2_index   = Channel.empty()
    ch_bowtie2_version = Channel.empty() 
    if ('bowtie2' in prepare_tool_indices) {
        if (params.host_bowtie2_index) {
            if (params.host_bowtie2_index.endsWith('.tar.gz')) {
                ch_host_bowtie2_index = UNTAR_BOWTIE2_INDEX ( params.host_bowtie2_index ).untar
            } else {
                ch_host_bowtie2_index = file(params.host_bowtie2_index)
            }
        } else if (params.host_fasta) {        
            ch_host_bowtie2_index   = BOWTIE2_BUILD ( ch_host_fasta ).index
            ch_bowtie2_version = BOWTIE2_BUILD.out.version
        }
    }

    /* Uncompress BWA index or generate from scratch if required */ 
    ch_host_bwa_index   = Channel.empty()
    ch_bwa_version = Channel.empty()
    if ('bwa' in prepare_tool_indices) {
        if (params.host_bwa_index) {
            if (params.host_bwa_index.endsWith('.tar.gz')) {
                ch_host_bwa_index = UNTAR_BWA_INDEX (params.host_bwa_index ).untar
            } else {
                ch_host_bwa_index = file(params.host_bwa_index)
            }
        } else if (params.host_fasta) {
            ch_host_bwa_index   = BWA_INDEX ( ch_host_fasta ).index
            ch_bwa_version = BWA_INDEX.out.version
        }
    }

//     
//     ch_host_bwa_index   = Channel.empty()
//     ch_bwa_version = Channel.empty()

//     /* Uncompress BOWTIE2 index or generate from scratch if required */  
//     ch_host_bowtie2_index   = Channel.empty()
//     ch_bowtie2_version = Channel.empty() 
    
//     if (!params.host_fasta) {
//         ch_host_fasta = Channel.empty()
//     } else if (params.host_fasta.endsWith('.gz')) {
//         ch_host_fasta = GUNZIP_FASTA ( params.host_fasta ).gunzip
//     } else {
//         ch_host_fasta = file(params.host_fasta)
//     } else if ('bwa' in prepare_tool_indices){
//         if (params.host_bwa_index) {
//             if (params.host_bwa_index.endsWith('.tar.gz')) {
//                 ch_host_bwa_index = UNTAR_BWA_INDEX ( params.host_bwa_index ).untar
//             } else {
//                 ch_host_bwa_index = file(params.host_bwa_index)
//             }
//         } else {
//             ch_host_bwa_index   = BWA_INDEX ( ch_host_fasta ).index
//             ch_bwa_version = BWA_INDEX.out.version
//         }
    // } else if ('bowtie2' in prepare_tool_indices) {
    //     if (params.host_bowtie2_index) {
    //         if (params.host_bowtie2_index.endsWith('.tar.gz')) {
    //             ch_host_bowtie2_index = UNTAR_BOWTIE2_INDEX ( params.host_bowtie2_index ).untar
    //         } else {
    //             ch_host_bowtie2_index = file(params.host_bowtie2_index)
    //         }
    //     } else {        
    //         ch_host_bowtie2_index   = BOWTIE2_BUILD ( ch_host_fasta ).index
    //         ch_bowtie2_version = BOWTIE2_BUILD.out.version
    //     }
    // }
    
    
    /* Uncompress BWA index or generate from scratch if required */  
    // ch_host_bwa_index   = Channel.empty()
    // ch_bwa_version = Channel.empty()
    // if ('bwa' in prepare_tool_indices) {
    //     if (params.host_bwa_index) {
    //         if (params.host_bwa_index.endsWith('.tar.gz')) {
    //             ch_host_bwa_index = UNTAR_BWA_INDEX ( params.host_bwa_index ).untar
    //         } else {
    //             ch_host_bwa_index = file(params.host_bwa_index)
    //         }
    //     } else {        
    //         ch_host_bwa_index   = BWA_INDEX ( ch_host_fasta ).index
    //         ch_bwa_version = BWA_INDEX.out.version
    //     }
    // }

    /* Uncompress BOWTIE2 index or generate from scratch if required */  
    // ch_host_bowtie2_index   = Channel.empty()
    // ch_bowtie2_version = Channel.empty() 
    // if ('bowtie2' in prepare_tool_indices) {
    //     if (params.host_bowtie2_index) {
    //         if (params.host_bowtie2_index.endsWith('.tar.gz')) {
    //             ch_host_bowtie2_index = UNTAR_BOWTIE2_INDEX ( params.host_bowtie2_index ).untar
    //         } else {
    //             ch_host_bowtie2_index = file(params.host_bowtie2_index)
    //         }
    //     } else {        
    //         ch_host_bowtie2_index   = BOWTIE2_BUILD ( ch_host_fasta ).index
    //         ch_bowtie2_version = BOWTIE2_BUILD.out.version
    //     }
    // }

    emit:
    fasta            = ch_host_fasta         // path: genome.fasta
    bwa_index        = ch_host_bwa_index     // path: bwa/index/
    bowtie2_index    = ch_host_bowtie2_index // path: bowtie2/index/
}