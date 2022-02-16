// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process SNPEFF_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::snpeff=5.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1'
    } else {
        container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'
    }

    input:
    tuple val(meta), path(vcf)
    path db
    path config
    path fasta

    output:
    tuple val(meta), path("*.snpEff.csv")   , emit: csv
    tuple val(meta), path("*.vcf.gz*")      , emit: vcf
    tuple val(meta), path("*.txt")          , emit: txt
    tuple val(meta), path("*.html")         , emit: html

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def lastPath   = params.fasta.lastIndexOf(File.separator)
    def lastExt    = params.fasta.lastIndexOf(".")
    def index_base = params.fasta.substring(lastPath+1,lastExt)

    """
    snpEff ${index_base} \\
        -config $config \\
        -dataDir $db \\
        $vcf \\
        -csvStats ${prefix}.snpEff.csv \\
        | bgzip -c > ${prefix}.snpEff.vcf.gz

    tabix -p vcf -f ${prefix}.snpEff.vcf.gz
    mv snpEff_summary.html ${prefix}.snpEff.summary.html

    SnpSift extractFields -s "," \\
        -e "." \\
        ${prefix}.snpEff.vcf.gz \\
        CHROM POS REF ALT \\
        "ANN[*].GENE" "ANN[*].GENEID" \\
        "ANN[*].IMPACT" "ANN[*].EFFECT" \\
        "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
        "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
        "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
        "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
        "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
        "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
        > ${prefix}.snpSift.table.txt
    	"""
}