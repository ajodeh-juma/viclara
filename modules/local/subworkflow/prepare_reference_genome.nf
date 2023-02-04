/*
 * Uncompress and prepare reference genome files
*/

params.genome_options         = [:]
params.index_options          = [:]
params.bwa_index_options      = [:]
params.bowtie2_index_options  = [:]


include { GUNZIP as GUNZIP_FASTA           }  from '../process/gunzip'                                     addParams( options: params.genome_options               )
include { UNTAR as UNTAR_BWA_INDEX         }  from '../process/untar'                                      addParams( options: params.bwa_index_options            )
include { UNTAR as UNTAR_BOWTIE2_INDEX     }  from '../process/untar'                                      addParams( options: params.bowtie2_index_options        )
include { BWA_INDEX                        }  from '../../nf-core/software/bwa/index/main'                 addParams( options: params.bwa_index_options            )
include { BOWTIE2_BUILD                    }  from '../../nf-core/software/bowtie2/build/main'             addParams( options: params.bowtie2_index_options        )

workflow PREPARE_REFERENCE_GENOME {
    take:
    prepare_tool_indices // list: tools to prepare indices for

    main:
    /*
     * Uncompress genome fasta file if required
     */
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    /*
     * Uncompress BWA index or generate from scratch if required
     */  
    ch_bwa_index   = Channel.empty()
    ch_bwa_version = Channel.empty()  
    if ('bwa' in prepare_tool_indices) {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( params.bwa_index ).untar
            } else {
                ch_bwa_index = file(params.bwa_index)
            }
        } else {        
            ch_bwa_index   = BWA_INDEX ( ch_fasta ).index
            ch_bwa_version = BWA_INDEX.out.version
        }
    }

    /*
     * Uncompress BOWTIE2 index or generate from scratch if required
     */  
    ch_bowtie2_index   = Channel.empty()
    ch_bowtie2_version = Channel.empty()  
    if ('bowtie2' in prepare_tool_indices) {
        if (params.bowtie2_index) {
            if (params.bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( params.bowtie2_index ).untar
            } else {
                ch_bowtie2_index = file(params.bowtie2_index)
            }
        } else {        
            ch_bowtie2_index   = BOWTIE2_BUILD ( ch_fasta ).index
            ch_bowtie2_version = BOWTIE2_BUILD.out.version
        }
    }

    emit:
    fasta            = ch_fasta         // path: genome.fasta
    bwa_index        = ch_bwa_index     // path: bwa/index/
    bowtie2_index    = ch_bowtie2_index // path: bowtie2/index/
}