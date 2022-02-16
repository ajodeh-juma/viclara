// Uncompress and prepare reference genome files


params.bwa_index_options      = [:]
params.bowtie2_index_options  = [:]

include { UNTAR         }  from '../../local/process/untar'                             addParams( options: [:]                          )
include { BWA_INDEX     }  from '../../nf-core/software/bwa/index/main'                 addParams( options:  params.bwa_index_options    )
include { BOWTIE2_BUILD }  from '../../nf-core/software/bowtie2/build/main'             addParams( options: params.bowtie2_index_options )

workflow PREPARE_PHIX_GENOME {
    take:
    prepare_tool_indices // list: tools to prepare indices for

    main:
    phix_fasta = 'RTA'? params.genomes[ 'RTA' ].fasta ?: false : false
    ch_phix_fasta = Channel
        .fromPath(phix_fasta, checkIfExists: true)
        .ifEmpty { exit 1, "PhiX Genome Fasta file not found: ${phix_fasta}, you can download the files from (http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz) and extract in the igenomes dir" }

    // Uncompress BWA index or generate from scratch if required

    ch_bwa_index   = Channel.empty()
    ch_bwa_version = Channel.empty()  
    if ('bwa' in prepare_tool_indices) {
        ch_bwa_index   = BWA_INDEX ( ch_phix_fasta ).index
        ch_bwa_version = BWA_INDEX.out.version
    }

    // Uncompress BOWTIE2 index or generate from scratch if required
       
    ch_bowtie2_index   = Channel.empty()
    ch_bowtie2_version = Channel.empty()  
    if ('bowtie2' in prepare_tool_indices) {
        ch_bowtie2_index   = BOWTIE2_BUILD ( ch_phix_fasta ).index
        ch_bowtie2_version = BOWTIE2_BUILD.out.version
    }

    emit:
    fasta            = ch_phix_fasta    // path: genome.fasta
    bwa_index        = ch_bwa_index     // path: bwa/index/
    bowtie2_index    = ch_bowtie2_index // path: bowtie2/index/
}