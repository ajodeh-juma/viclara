// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNPEFF_BUILD_DB {
    tag "$fasta"
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
    path fasta
    path gff

    output:
    path("snpeff_db"), emit: db
    path("*.config"),  emit: config
        
    script:
    def software = getSoftwareName(task.process)
    def lastPath   = fasta.lastIndexOf(File.separator)
    def lastExt    = fasta.lastIndexOf(".")
    def index_base = params.fasta.substring(lastPath+1,lastExt)
        
    """
    echo "${index_base}.genome : ${index_base}" > snpeff.config   
    build_snpeff_db.py --reference $fasta --gff $gff --config snpeff.config --snpeff-db snpeff_db
    echo \$(snpEff -version 2>&1) | sed 's/^.*SnpEff //; s/ .*\$//' > ${software}.version.txt
    """
    
}