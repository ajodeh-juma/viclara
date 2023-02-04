////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////
nextflow.enable.dsl = 2

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.metadata, params.multiqc_config,
    params.host_fasta, params.host_bwa_index, params.host_bowtie2_index
]
for (param in checkPathParamList) { 
    if (param) { 
        file(param, checkIfExists: true) 
    } 
}

def segmentsList = ['S', 'M', 'L']
if (!segmentsList.contains(params.segment)) {
    exit 1, "Invalid segment option: ${params.segment}. Valid options: ${segmentsList.join(', ')}"
} else {
    ch_prefix = params.segment + '-Segment'
}

// Check metadata file
if (!params.metadata) {
    ch_metadata = file("$projectDir/assets/metadata.csv", checkIfExists: true)
} else {
    ch_metadata = file(params.metadata)
}


// Check trimming tool
def trimmersList = ['fastp', 'trimmomatic']
if (!params.skip_trimming) {
    if (!trimmersList.contains(params.trimmer)) {
        exit 1, "Invalid trimmer option: ${params.trimmer}. Valid options: ${trimmersList.join(', ')}"
    }
}

// Check alignment tools
def prepareToolIndices  = []
def alignerList         = ['bwa', 'bowtie2']
if (!params.skip_alignment) {
    if (!alignerList.contains(params.aligner)) {
        exit 1, "Invalid aligner option: ${params.aligner}. Valid options: ${alignerList.join(', ')}"
    }
    prepareToolIndices << params.aligner
}

// Check variant calling tools
def variant_callersList = ['bcftools', 'varscan']
if (!variant_callersList.contains(params.variant_caller)) {
    exit 1, "Invalid variant calling option: ${params.variant_caller}. Valid options: ${variant_callersList.join(', ')}"
}


// check database
if (!params.skip_kraken2 && params.database) {
    ch_database = Channel.fromPath(params.database, type: 'dir')
} else if (params.skip_kraken2 && params.database == null) {
    exit 1, "Missing options, database path: ${params.database}"
}


////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]


// add options to multiqc module
def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// add options to bwa index module
def bwa_index_options     = modules['bwa_index']
if (!params.save_reference)  { bwa_index_options['publish_files'] = false }

// add options to bwa mem module
def bwa_mem_options         = modules['bwa_mem']
if (params.save_align_intermeds) { bwa_mem_options.publish_files.put('bam','') }

// add options to bowtie2 build module
def bowtie2_build_options     = modules['bowtie2_build']
if (!params.save_reference)  { bowtie2_build_options['publish_files'] = false }

// add options to bowtie2 align module
def bowtie2_align_options         = modules['bowtie2_align']
if (params.save_align_intermeds) { bowtie2_align_options.publish_files.put('bam','') }

def samtools_sort_options_phix = modules['samtools_sort_phix']
if (['bwa','bowtie2'].contains(params.aligner)) {
    if (params.save_align_intermeds && params.skip_markduplicates) {
        samtools_sort_options.publish_files.put('bam','')
        samtools_sort_options.publish_files.put('bai','')
    }
}

def samtools_sort_options_host = modules['samtools_sort_host']
if (['bwa','bowtie2'].contains(params.aligner)) {
    if (params.save_align_intermeds && params.skip_markduplicates) {
        samtools_sort_options.publish_files.put('bam','')
        samtools_sort_options.publish_files.put('bai','')
    }
}

def samtools_sort_options_reference = modules['samtools_sort_reference']
if (['bwa','bowtie2'].contains(params.aligner)) {
    if (params.save_align_intermeds && params.skip_markduplicates) {
        samtools_sort_options.publish_files.put('bam','')
        samtools_sort_options.publish_files.put('bai','')
    }
}


// Get total number of mapped reads from flagstat file
def get_mapped_from_flagstat(flagstat) {
    def mapped = 0
    flagstat.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped = line.tokenize().first().toInteger()
        }
    }
    return mapped
}

// Function that checks the number of mapped reads from flagstat output and returns true if > params.min_mapped_reads and otherwise false
pass_mapped_reads = [:]
fail_mapped_reads = [:]
def check_mapped(sample, flagstat, min_mapped=500) {
    mapped = get_mapped_from_flagstat(flagstat)
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    if (mapped < min_mapped.toInteger()) {
        log.info ">${c_red}>>>> $sample FAILED MAPPED READ THRESHOLD: ${mapped} < ${params.min_mapped}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! <<<<${c_reset}<"
        fail_mapped_reads[sample] = mapped
        return false
    } else {
        pass_mapped_reads[sample] = mapped
        return true
    }
}

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { MULTIQC                                     } from './modules/local/process/multiqc'                               addParams( options: multiqc_options                                                                                                                                                           )
include { GET_SOFTWARE_VERSIONS                       } from './modules/local/process/get_software_versions'                 addParams( options: [publish_files : ['csv':'']]                                                                                                                                              )
include { TRIMMOMATIC                                 } from './modules/local/process/trimmomatic'                           addParams( options: [:]                                                                                                                                                                       )
include { PREPARE_PHIX_GENOME                         } from './modules/local/subworkflow/prepare_phix_genome'               addParams( bwa_index_options: bwa_index_options, bowtie2_build_options: bowtie2_build_options                                                                                                 )
include { PREPARE_GENOME                              } from './modules/local/subworkflow/prepare_genome'                    addParams( bwa_index_options: bwa_index_options, bowtie2_build_options: bowtie2_build_options                                                                                                 )
include { PREPARE_HOST_GENOME                         } from './modules/local/subworkflow/prepare_host_genome'               addParams( genome_options: publish_genome_options, index_options: publish_index_options, bwa_index_options: bwa_index_options, bowtie2_index_options: bowtie2_build_options                   )
include { ALIGN_BWA as ALIGN_BWA_PHIX                 } from './modules/local/subworkflow/align_bwa'                         addParams( bwa_mem_options: bwa_mem_options, samtools_options: samtools_sort_options_phix                                                                                                     )
include { ALIGN_BWA as ALIGN_BWA_HOST                 } from './modules/local/subworkflow/align_bwa'                         addParams( bwa_mem_options: bwa_mem_options, samtools_options: samtools_sort_options_host                                                                                                     )
include { ALIGN_BWA as ALIGN_BWA_FASTA                } from './modules/local/subworkflow/align_bwa'                         addParams( align_options: bwa_mem_options, samtools_options: samtools_sort_options_reference                                                                                                  )
include { ALIGN_BOWTIE2 as ALIGN_BOWTIE2_PHIX         } from './modules/local/subworkflow/align_bowtie2'                     addParams( align_options: bowtie2_align_options, samtools_options: samtools_sort_options_phix                                                                                                 )
include { ALIGN_BOWTIE2 as ALIGN_BOWTIE2_HOST         } from './modules/local/subworkflow/align_bowtie2'                     addParams( align_options: bowtie2_align_options, samtools_options: samtools_sort_options_host                                                                                                 )
include { ALIGN_BOWTIE2 as ALIGN_BOWTIE2_FASTA        } from './modules/local/subworkflow/align_bowtie2'                     addParams( align_options: bowtie2_align_options, samtools_options: samtools_sort_options_reference                                                                                            )
include { SAMTOOLS_UNMAPPED_FASTQ as PHIX_UNMAPPED    } from './modules/local/process/samtools_unmapped_fastq'               addParams( options: modules['samtools_unmapped_fastq_phix']                                                                                                                                   )
include { SAMTOOLS_UNMAPPED_FASTQ as HOST_UNMAPPED    } from './modules/local/process/samtools_unmapped_fastq'               addParams( options: modules['samtools_unmapped_fastq_host']                                                                                                                                   )
include { SAMTOOLS_MAPPED_BAM                         } from './modules/local/process/samtools_mapped_bam'                   addParams( options: modules['samtools_mapped_bam']                                                                                                                                            )
include { MARKDUPLICATES_PICARD                       } from './modules/local/subworkflow/markduplicates_picard'             addParams( picard_markduplicates_options: modules['picard_markduplicates'], samtools_options: modules['picard_markduplicates_samtools']                                                       )
include { BEDTOOLS_GENOME_COVERAGE                    } from './modules/local/process/bedtools_genome_coverage'              addParams( options: modules['bedtools_genomecov']                                                                                                                                             )
include { COVERAGE_PLOTS                              } from './modules/local/process/coverage_plots'                        addParams( options: modules['plots']                                                                                                                                                 )
include { VARSCAN_MPILEUP                             } from './modules/local/process/varscan_mpileup'                       addParams(options: modules['varscan_mpileup']                                                                                                                                                 )
include { VARSCAN_CONSENSUS                           } from './modules/local/process/varscan_consensus'                     addParams(options: modules['varscan_consensus']                                                                                                                                               )
include { BCFTOOLS_CONSENSUS                          } from './modules/local/process/bcftools_consensus'                    addParams(options: modules['bcftools_consensus']                                                                                                                                              )
include { SNPEFF_BUILD_DB                             } from './modules/local/process/snpeff_build_db'                       addParams(options: [:]                                                                                                                                                                        )
include { SNPEFF_ANNOTATE as SNPEFF_VARSCAN_LOFREQ
          SNPEFF_ANNOTATE as SNPEFF_VARSCAN_HIFREQ
          SNPEFF_ANNOTATE as SNPEFF_BCFTOOLS          } from './modules/local/process/snpeff_annotate'                       addParams(options: modules['snpeff_annotate']                                                                                                                                                 )
include { TRINITY_ASSEMBLE                            } from './modules/local/process/trinity_assemble'                      addParams(options: [:])
include { METASPADES_ASSEMBLE                         } from './modules/local/process/metaspades_assemble'                   addParams(options: modules['metaspades_assemble']                                                                                                                                             )
include { TADPOLE_ASSEMBLE                            } from './modules/local/process/tadpole_assemble'                      addParams(options: [:]                                                                                                                                                                        )
include { ABACAS_CONTIGUATE                           } from './modules/local/process/abacas_contiguate'                     addParams(options: modules['abacas_contiguate']                                                                                                                                               )
include { KRAKEN2_KRAKEN2                             } from './modules/local/process/kraken2.nf'                            addParams(options: [:]                                                                                                                                                                        )
include { KRAKEN2_MPA_SUMMARY                         } from './modules/local/process/kraken2_mpa_summary.nf'                addParams(options: modules['summary_stats']                                                                                                                                                                        )
include { KRAKEN2_MPA_PLOT                            } from './modules/local/process/kraken2_mpa_plot.nf'                   addParams(options: modules['plots']                                                                                                                                                                        )
include { CONSENSUS_QC                                } from './modules/local/process/consensus_qc'                          addParams( options: modules['consensus_qc']                                                                                                                                                   )
include { QC_SUMMARY_CSV                              } from './modules/local/process/qc_summary_csv'                        addParams( options: modules['consensus_qc']                                                                                                                                                   )
include { READ_COUNT as READS_RAW                     } from './modules/local/process/read_count'                            addParams( options: modules['raw_counts']                                                                                                                                                      )
include { READ_COUNT as READS_TRIMMED                 } from './modules/local/process/read_count'                            addParams( options: modules['trimmed_counts']                                                                                                                                                  )
include { PARSE_FLAGSTAT                              } from './modules/local/process/parse_alignment_stats'                 addParams( options: modules['alignment_stats']                                                                                                                                               )
include { MERGE_COUNTS as MERGE_RAW_COUNTS            } from './modules/local/process/merge_counts'                          addParams( options: modules['merge_raw_counts']                                                                                                                                                )
include { MERGE_COUNTS as MERGE_TRIMMED_COUNTS        } from './modules/local/process/merge_counts'                          addParams( options: modules['merge_trimmed_counts']                                                                                                                                                 )
include { MERGE_COUNTS as MERGE_MAPPED_COUNTS         } from './modules/local/process/merge_counts'                          addParams( options: modules['merge_mapped_counts']                                                                                                                                               )
include { MERGE_STATS                                 } from './modules/local/process/merge_stats'                           addParams( options: modules['summary_stats']                                                                                                                                               )
include { PLOT_READ_COUNTS                            } from './modules/local/process/plot_read_counts'                      addParams( options: modules['plots']                                                                                                                                               )


////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

// MODULE: Installed directly from nf-core/modules

include { FASTQC                     } from './modules/nf-core/software/fastqc/main'                 addParams( options: modules['fastqc']           )
include { FASTP                      } from './modules/nf-core/software/fastp/main'                  addParams( options: modules['fastp']            )
include { BCFTOOLS_MPILEUP           } from './modules/nf-core/software/bcftools/mpileup/main'       addParams( options: modules['bcftools_mpileup'] )
include { SAMTOOLS_MPILEUP           } from './modules/nf-core/software/samtools/mpileup/main'       addParams( options: modules['samtools_mpileup'] )



////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []


// Create a channel for input read files
if (params.input_paths) {
    if (params.single_end) {
        ch_reads = Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [file(row[1][0], checkIfExists: true)]]}
            .dump()
            .ifEmpty {exit 1, "$params.input_paths was empty - no input files supplied"}
    } else {
        ch_reads = Channel
            .from(params.input_paths)
            .map {row -> [row[0], [file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true)]]}
            .dump()
            .ifEmpty {exit 1, "$params.input_paths was empty - no input files supplied"}
    }
} else {
    if (params.single_end) {
        ch_reads = Channel
            .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
            .map { row -> 
            def meta = [:]
            meta.id = row[0]
            [ meta , [ file(row[1][0], checkIfExists: true) ] ]}
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --single_end on the command line." }
    } else {
        ch_reads = Channel
            .fromFilePairs(params.input, size: params.single_end ? 1 : 2 )
            .map { row -> 
            def meta = [:]
            meta.id = row[0]
            [ meta , [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ]}
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --single_end on the command line." }
            }
}



workflow viclara {
    
    // channel to collect software versions
    ch_software_versions = Channel.empty()

    // println ch_reads.view()

    // MODULE: FastQC
    ch_fastqc_multiqc = Channel.empty()
    if (!params.skip_qc) {
        FASTQC (ch_reads)
        ch_fastqc_multiqc = FASTQC.out.zip
        ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

        // count raw reads
        // ch_raw_counts = RAW_COUNTS(ch_reads)

        READS_RAW(ch_reads)
        ch_raw_counts = READS_RAW.out.csv
    }

    // MODULE: Trim reads
    if (params.skip_trimming) {
        ch_trimmed_reads = ch_reads
        ch_trimmed_counts = ch_raw_counts
    }

    ch_fastp_multiqc = Channel.empty()
    if (!params.skip_trimming && params.trimmer == 'fastp') {
        FASTP (ch_reads)
        ch_trimmed_reads     = FASTP.out.reads
        ch_fastp_multiqc     = FASTP.out.log
        ch_software_versions = ch_software_versions.mix(FASTP.out.version.first().ifEmpty(null))

        // count trimmed reads
        // ch_trimmed_counts = TRIMMED_COUNTS(ch_trimmed_reads)
        READS_TRIMMED(ch_trimmed_reads)
        ch_trimmed_counts    = READS_TRIMMED.out.csv
    }
    
    ch_trimmomatic_multiqc = Channel.empty()
    if (!params.skip_trimming && params.trimmer == 'trimmomatic' ) {
        TRIMMOMATIC(ch_reads)
        ch_trimmed_reads     = TRIMMOMATIC.out.reads
        ch_trimmomatic_multiqc = TRIMMOMATIC.out.log
        ch_software_versions = ch_software_versions.mix(TRIMMOMATIC.out.version.first().ifEmpty(null))

        // count trimmed reads
        // ch_trimmed_counts = TRIMMED_COUNTS(ch_trimmed_reads)
        READS_TRIMMED(ch_trimmed_reads)
        ch_trimmed_counts    = READS_TRIMMED.out.csv
    }

    // SUBWORKFLOW: classify reads with Kraken2 and visualize
    if (!params.skip_kraken2 && params.database) {
        KRAKEN2_KRAKEN2(ch_trimmed_reads, ch_database.collect())

        ch_reports_txt = KRAKEN2_KRAKEN2.out.txt
        ch_reports_txt
            .filter { row -> file(row[1])}
            .map { it[1] }
            .set { ch_reports_txt }
        KRAKEN2_MPA_SUMMARY(ch_reports_txt.collect(), ch_prefix)
        KRAKEN2_MPA_PLOT(KRAKEN2_MPA_SUMMARY.out.species, ch_prefix)
    }

    // SUBWORKFLOW: Filter PhiX reads and extract unmapped reads
    if (!params.filter_phix) {
        ch_bwa_phix_multiqc     = Channel.empty()
        ch_bowtie2_phix_multiqc = Channel.empty()
        ch_phix_filtered_reads  = ch_trimmed_reads
    } else if (params.filter_phix && params.aligner == 'bwa') {
       PREPARE_PHIX_GENOME(prepareToolIndices)
        ALIGN_BWA_PHIX (
            ch_trimmed_reads,
            PREPARE_PHIX_GENOME.out.bwa_index.collect()
        )
        ch_genome_bam        = ALIGN_BWA_PHIX.out.bam
        ch_genome_bai        = ALIGN_BWA_PHIX.out.bai
        ch_samtools_stats    = ALIGN_BWA_PHIX.out.stats
        ch_samtools_flagstat = ALIGN_BWA_PHIX.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA_PHIX.out.idxstats
        PHIX_UNMAPPED(ch_genome_bam)
        ch_phix_filtered_reads  = PHIX_UNMAPPED.out.reads
    } else if (params.filter_phix && params.aligner == 'bowtie2') {
        PREPARE_PHIX_GENOME(prepareToolIndices)
        ALIGN_BOWTIE2_PHIX (
            ch_trimmed_reads,
            PREPARE_PHIX_GENOME.out.bowtie2_index.collect()
        )
        ch_genome_bam        = ALIGN_BOWTIE2_PHIX.out.bam
        ch_genome_bai        = ALIGN_BOWTIE2_PHIX.out.bai
        ch_bowtie2_multiqc   = ALIGN_BOWTIE2_PHIX.out.log
        ch_samtools_stats    = ALIGN_BOWTIE2_PHIX.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2_PHIX.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2_PHIX.out.idxstats
        PHIX_UNMAPPED(ch_genome_bam)
        ch_phix_filtered_reads  = PHIX_UNMAPPED.out.reads
    }

    // SUBWORKFLOW: Filter host reads and extract unmapped reads
    if (!params.host_fasta || !params.host_bwa_index || !params.host_bowtie2_index) {
        ch_bwa_phix_multiqc     = Channel.empty()
        ch_bowtie2_phix_multiqc = Channel.empty()
        ch_host_filtered_reads  = ch_phix_filtered_reads
    }
    if ( params.host_fasta && params.aligner == 'bwa') {
        PREPARE_HOST_GENOME(prepareToolIndices)
        ALIGN_BWA_HOST (
            ch_phix_filtered_reads,
            PREPARE_HOST_GENOME.out.bwa_index
        )
        ch_genome_bam        = ALIGN_BWA_HOST.out.bam
        ch_genome_bai        = ALIGN_BWA_HOST.out.bai
        ch_samtools_stats    = ALIGN_BWA_HOST.out.stats
        ch_samtools_flagstat = ALIGN_BWA_HOST.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA_HOST.out.idxstats
        HOST_UNMAPPED(ch_genome_bam)
        ch_host_filtered_reads  = HOST_UNMAPPED.out.reads
    } 
    if ( params.host_bwa_index ) {
        PREPARE_HOST_GENOME(prepareToolIndices)
        ALIGN_BWA_HOST (
            ch_phix_filtered_reads,
            PREPARE_HOST_GENOME.out.bwa_index
        )
        ch_genome_bam        = ALIGN_BWA_HOST.out.bam
        ch_genome_bai        = ALIGN_BWA_HOST.out.bai
        ch_samtools_stats    = ALIGN_BWA_HOST.out.stats
        ch_samtools_flagstat = ALIGN_BWA_HOST.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA_HOST.out.idxstats
        HOST_UNMAPPED(ch_genome_bam)
        ch_host_filtered_reads  = HOST_UNMAPPED.out.reads
    }
    if (params.host_fasta && params.aligner == 'bowtie2') {
        PREPARE_HOST_GENOME(prepareToolIndices)
        ALIGN_BOWTIE2_HOST (
            ch_phix_filtered_reads,
            PREPARE_HOST_GENOME.out.bowtie2_index.collect()
        )
        ch_genome_bam        = ALIGN_BOWTIE2_HOST.out.bam
        ch_genome_bai        = ALIGN_BOWTIE2_HOST.out.bai
        ch_bowtie2_multiqc   = ALIGN_BOWTIE2_HOST.out.log
        ch_samtools_stats    = ALIGN_BOWTIE2_HOST.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2_HOST.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2_HOST.out.idxstats
        HOST_UNMAPPED(ch_genome_bam)
        ch_host_filtered_reads  = HOST_UNMAPPED.out.reads
    }
    if ( params.host_bowtie2_index ) {
        PREPARE_HOST_GENOME(prepareToolIndices)
        ALIGN_BOWTIE2_HOST (
            ch_phix_filtered_reads,
            PREPARE_HOST_GENOME.out.bowtie2_index
        )
        ch_genome_bam        = ALIGN_BOWTIE2_HOST.out.bam
        ch_genome_bai        = ALIGN_BOWTIE2_HOST.out.bai
        ch_bowtie2_multiqc   = ALIGN_BOWTIE2_HOST.out.log
        ch_samtools_stats    = ALIGN_BOWTIE2_HOST.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2_HOST.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2_HOST.out.idxstats
        HOST_UNMAPPED(ch_genome_bam)
        ch_host_filtered_reads  = HOST_UNMAPPED.out.reads
    }

    // SUBWORKFLOW: Alignment to the reference FASTA genome with BWA MEM
    ch_bwa_multiqc = Channel.empty()
    ch_sort_bam    = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'bwa') {
        PREPARE_GENOME(prepareToolIndices)
        ALIGN_BWA_FASTA (
            ch_host_filtered_reads,
            PREPARE_GENOME.out.bwa_index
        )
        ch_fasta             = PREPARE_GENOME.out.fasta
        ch_gff               = PREPARE_GENOME.out.gff
        ch_genome_bam        = ALIGN_BWA_FASTA.out.bam
        ch_genome_bai        = ALIGN_BWA_FASTA.out.bai
        ch_samtools_stats    = ALIGN_BWA_FASTA.out.stats
        ch_samtools_flagstat = ALIGN_BWA_FASTA.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA_FASTA.out.idxstats
    }

    // SUBWORKFLOW: Alignment to the reference FASTA genome with BOWTIE2
    ch_bowtie2_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'bowtie2') {
        PREPARE_GENOME(prepareToolIndices)    
        ALIGN_BOWTIE2_FASTA (
            ch_host_filtered_reads,
            PREPARE_GENOME.out.bowtie2_index
        )
        ch_fasta             = PREPARE_GENOME.out.fasta
        ch_gff               = PREPARE_GENOME.out.gff
        ch_genome_bam        = ALIGN_BOWTIE2_FASTA.out.bam
        ch_genome_bai        = ALIGN_BOWTIE2_FASTA.out.bai
        ch_bowtie2_multiqc   = ALIGN_BOWTIE2_FASTA.out.log
        ch_samtools_stats    = ALIGN_BOWTIE2_FASTA.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2_FASTA.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2_FASTA.out.idxstats
        
    }

    // add flagstat channel 
    ch_sort_bam = ch_genome_bam.join(ch_samtools_flagstat)
    ch_alignment_stats = PARSE_FLAGSTAT(ch_samtools_flagstat)

    // MODULE: merge the read counts and mapping stats
    ch_raw_counts
        .filter { row -> file(row[1]) }
        .map { it[1] }
        .set { ch_raw_counts }

    ch_trimmed_counts
        .filter { row -> file(row[1]) }
        .map { it[1] }
        .set { ch_trimmed_counts }

    ch_alignment_stats
        .filter { row -> file(row[1]) }
        .map { it[1] }
        .set { ch_alignment_stats }

    MERGE_RAW_COUNTS ( ch_raw_counts.collect(), ch_prefix )
    MERGE_TRIMMED_COUNTS ( ch_trimmed_counts.collect(), ch_prefix )
    MERGE_MAPPED_COUNTS ( ch_alignment_stats.collect(), ch_prefix )

    ch_raw_rc = MERGE_RAW_COUNTS.out.csv
    ch_trimmed_rc = MERGE_TRIMMED_COUNTS.out.csv
    ch_mapped_rc  = MERGE_MAPPED_COUNTS.out.csv

    MERGE_STATS (ch_raw_rc, ch_trimmed_rc, ch_mapped_rc, ch_prefix)
    PLOT_READ_COUNTS (MERGE_STATS.out.csv, ch_prefix)
    

    // Remove samples that failed mapped read threshold
    ch_sort_bam
        .filter { row -> check_mapped( row[0], row[2], params.min_mapped )  }
        .map {it[0..1]}
        .set { ch_sort_bam }

    // SUBWORKFLOW: Extract mapped alignments and index
    SAMTOOLS_MAPPED_BAM (ch_sort_bam)

    // SUBWORKFLOW: Mark duplicates 
    if (!params.skip_markduplicates) {
        MARKDUPLICATES_PICARD(SAMTOOLS_MAPPED_BAM.out.bam)
        ch_bam = MARKDUPLICATES_PICARD.out.bam
    } else {
        ch_bam = ch_sort_bam
    }

    // MODULE: Genome wide coverage
    BEDTOOLS_GENOME_COVERAGE (ch_bam)
    
    // MODULE: Call variants, Generate consensus sequence
    if (params.variant_caller == 'bcftools') {
        BCFTOOLS_MPILEUP(ch_bam, ch_fasta.collect())
        BCFTOOLS_CONSENSUS(BCFTOOLS_MPILEUP.out.vcf, BCFTOOLS_MPILEUP.out.tbi, ch_bam, ch_fasta.collect())
        ch_consensus_fasta = BCFTOOLS_CONSENSUS.out.fasta

        // MODULE: Build SnpEff db and Annotate variants 
        SNPEFF_BCFTOOLS(BCFTOOLS_MPILEUP.out.vcf, ch_fasta.collect(), ch_gff.collect())
    }
    if (params.variant_caller == 'varscan') {
        SAMTOOLS_MPILEUP(ch_bam, ch_fasta.collect())
        VARSCAN_MPILEUP(SAMTOOLS_MPILEUP.out.mpileup)

        ch_varscan_hifreq_vcf = VARSCAN_MPILEUP.out.vcf.map { row -> [row[0], [ file(row[1][1], checkIfExists: true) ] ] }
        ch_varscan_hifreq_tbi = VARSCAN_MPILEUP.out.tbi.map { row -> [row[0], [ file(row[1][1], checkIfExists: true) ] ] }

        ch_varscan_lofreq_vcf = VARSCAN_MPILEUP.out.vcf.map { row -> [row[0], [ file(row[1][0], checkIfExists: true) ] ] }
        ch_varscan_lofreq_tbi = VARSCAN_MPILEUP.out.tbi.map { row -> [row[0], [ file(row[1][0], checkIfExists: true) ] ] }

        VARSCAN_CONSENSUS(ch_varscan_hifreq_vcf, ch_varscan_hifreq_tbi, ch_bam, ch_fasta.collect())
        ch_consensus_fasta = VARSCAN_CONSENSUS.out.fasta

        // MODULE: Build SnpEff db and Annotate variants
        SNPEFF_VARSCAN_HIFREQ(ch_varscan_hifreq_vcf, ch_fasta.collect(), ch_gff.collect())
        SNPEFF_VARSCAN_LOFREQ(ch_varscan_lofreq_vcf, ch_fasta.collect(), ch_gff.collect())

    }

    // MODULE: generate consensus qc metrics

    // channel for all the consensus fasta files
    ch_consensus_fasta
        .map { row -> 
            def meta = [:]
            meta.id = row[0].id
            [meta , file(row[1][1])]
            }
        .set {ch_consensus_fa}

    // println (ch_consensus_fasta.view())
    
    ch_consensus_qc = ch_bam.join(ch_consensus_fa)
    CONSENSUS_QC (ch_consensus_qc, ch_fasta.collect())

    // MODULE: collate the qc csv files
    ch_qc_csv = CONSENSUS_QC.out.csv
    ch_qc_csv
        .filter { row -> file(row[1])}
        .map { it[1] }
        .set { ch_qc_csv }

    QC_SUMMARY_CSV (ch_qc_csv.collect(), MERGE_STATS.out.csv, ch_prefix)

    // MODULE: Assemble reads
    if (params.skip_assembly) {
        if (params.assembler == 'trinity') {
            TRINITY_ASSEMBLE(ch_host_filtered_reads)
            ABACAS_CONTIGUATE(TRINITY_ASSEMBLE.out.contig, ch_fasta.collect())
            ch_masked_fa = ABACAS_CONTIGUATE.out.fasta
        }
        if (params.assembler == 'metaspades') {
            METASPADES_ASSEMBLE(ch_host_filtered_reads)
            ABACAS_CONTIGUATE(METASPADES_ASSEMBLE.out.contigs, ch_fasta.collect())
            ch_masked_fa = ABACAS_CONTIGUATE.out.fasta
        }
        if (params.assembler == 'tadpole') {
            TADPOLE_ASSEMBLE(ch_host_filtered_reads)
            ABACAS_CONTIGUATE(TADPOLE_ASSEMBLE.out.contigs, ch_fasta.collect())
            ch_masked_fa = ABACAS_CONTIGUATE.out.fasta
        }
    }
    
    

    // channel for all the consensus fasta files
    ch_consensus_fasta
        .filter { row -> file(row[1][1])}
        .map { it[1][1] }
        .set { ch_consensus_fasta }

    // channel for all the genomewide coverage files
    ch_genome_coverage_bed = BEDTOOLS_GENOME_COVERAGE.out.coverage_bed
    ch_genome_coverage_bed
        .filter { row -> file(row[1])}
        .map { it[1] }
        .set { ch_genome_coverage_bed }

    COVERAGE_PLOTS (
        ch_consensus_fasta.collect(),
        ch_genome_coverage_bed.collect(),
        ch_metadata, 
        QC_SUMMARY_CSV.out.csv,
        ch_prefix
    )

    // MODULE: Pipeline reporting
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map { it }.collect())

    // MultiQC
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fastqc_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastp_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trimmomatic_multiqc.collect{it[1]}.ifEmpty([])
            //ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}


////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
