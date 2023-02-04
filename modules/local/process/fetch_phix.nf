// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process wget_PhiX {
    tag params.phix_url
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::wget=1.20.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gnu-wget1.18--0"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }
    
    output:
    path "*.tar.gz"                      , emit: phix_tar_gz
    path "*PhiX"                         , emit: base
    //path "*BWAIndex"                     , emit: bwaindex
    path "*.version.txt"                 , emit: version
    
    script:
    def software = getSoftwareName(task.process)
    def lastPath = "${params.phix_url}".lastIndexOf(File.separator)
    def base = "${params.phix_url}".substring(lastPath+1)
    // 'wget_PhiX' {
    //         publish_dir  = "${params.igenomes_base}"
    //     }


    """
    wget --continue ${params.phix_url} 
    tar -xzvf ${base}
    echo \$(wget --version 2>&1) | sed 's/^.*(GNU Wget) //; s/ Copyright.*\$//' > ${software}.version.txt
    """
}