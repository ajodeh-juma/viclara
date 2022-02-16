// Uncompress and prepare reference genome files


params.bwa_index_options      = [:]
params.bowtie2_index_options  = [:]

include { UNTAR         }  from '../../local/process/untar'                             addParams(options: [:])
include { BWA_INDEX     }  from '../../nf-core/software/bwa/index/main'                 addParams( options:  params.bwa_index_options     )
include { BOWTIE2_BUILD }  from '../../nf-core/software/bowtie2/build/main'             addParams( options: params.bowtie2_index_options    )

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list: tools to prepare indices for

    main:

    if (params.segment == 'S') {
        fasta = 'ZH-548'? params.genomes[ 'ZH-548' ][ 'S' ].fasta ?: false : false
        ch_fasta = Channel
            .fromPath(fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${fasta}, you can download the file from https://www.ncbi.nlm.nih.gov/nuccore/NC_014395.1?report=fasta"}

        gff = 'ZH-548'? params.genomes[ 'ZH-548' ][ 'S' ].gff ?: false : false
        ch_gff = Channel
            .fromPath(gff, checkIfExists: true)
            .ifEmpty { exit 1, "Genome annotation file not found: ${gff}"}

        
        // Uncompress BWA index or generate from scratch if required

        ch_bwa_index   = Channel.empty()
        ch_bwa_version = Channel.empty()  
        if ('bwa' in prepare_tool_indices) {
            ch_bwa_index   = BWA_INDEX ( fasta ).index
            ch_bwa_version = BWA_INDEX.out.version
        }

        // Uncompress BOWTIE2 index or generate from scratch if required
        
        ch_bowtie2_index   = Channel.empty()
        ch_bowtie2_version = Channel.empty()  
        if ('bowtie2' in prepare_tool_indices) {
            ch_bowtie2_index   = BOWTIE2_BUILD ( fasta ).index
            ch_bowtie2_version = BOWTIE2_BUILD.out.version
        }
    }

    if (params.segment == 'M') {
        fasta = 'ZH-548'? params.genomes[ 'ZH-548' ][ 'M' ].fasta ?: false : false
        ch_fasta = Channel
            .fromPath(fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${fasta}, you can download the file from https://www.ncbi.nlm.nih.gov/nuccore/NC_014396.1?report=fasta"}
        
        gff = 'ZH-548'? params.genomes[ 'ZH-548' ][ 'M' ].gff ?: false : false
        ch_gff = Channel
            .fromPath(gff, checkIfExists: true)
            .ifEmpty { exit 1, "Genome annotation file not found: ${gff}"}

        
        // Uncompress BWA index or generate from scratch if required

        ch_bwa_index   = Channel.empty()
        ch_bwa_version = Channel.empty()  
        if ('bwa' in prepare_tool_indices) {
            ch_bwa_index   = BWA_INDEX ( fasta ).index
            ch_bwa_version = BWA_INDEX.out.version
        }

        // Uncompress BOWTIE2 index or generate from scratch if required
        
        ch_bowtie2_index   = Channel.empty()
        ch_bowtie2_version = Channel.empty()  
        if ('bowtie2' in prepare_tool_indices) {
            ch_bowtie2_index   = BOWTIE2_BUILD ( fasta ).index
            ch_bowtie2_version = BOWTIE2_BUILD.out.version
        }
    }

    if (params.segment == 'L') {
        fasta = 'ZH-548'? params.genomes[ 'ZH-548' ][ 'L' ].fasta ?: false : false
        ch_fasta = Channel
            .fromPath(fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${fasta}, you can download the file from https://www.ncbi.nlm.nih.gov/nuccore/NC_014397.1?report=fasta"}

        gff = 'ZH-548'? params.genomes[ 'ZH-548' ][ 'L' ].gff ?: false : false
        ch_gff = Channel
            .fromPath(gff, checkIfExists: true)
            .ifEmpty { exit 1, "Genome annotation file not found: ${gff}"}

        // Uncompress BWA index or generate from scratch if required

        ch_bwa_index   = Channel.empty()
        ch_bwa_version = Channel.empty()  
        if ('bwa' in prepare_tool_indices) {
            ch_bwa_index   = BWA_INDEX ( fasta ).index
            ch_bwa_version = BWA_INDEX.out.version
        }

        // Uncompress BOWTIE2 index or generate from scratch if required
        
        ch_bowtie2_index   = Channel.empty()
        ch_bowtie2_version = Channel.empty()  
        if ('bowtie2' in prepare_tool_indices) {
            ch_bowtie2_index   = BOWTIE2_BUILD ( fasta ).index
            ch_bowtie2_version = BOWTIE2_BUILD.out.version
        }
    }
    emit:
    fasta            = ch_fasta         // path: genome.fasta
    gff              = ch_gff           // path: genome.gff
    bwa_index        = ch_bwa_index     // path: bwa/index/
    bowtie2_index    = ch_bowtie2_index // path: bowtie2/index/
}