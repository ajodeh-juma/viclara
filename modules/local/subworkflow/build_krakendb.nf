#!/usr/bin/env nextflow

/******************************************************************************
*************************** build kraken2 databases ***************************
******************************************************************************/

if (params.classify) {
    process build_krakendatabase {
        label 'high_memory'
        tag "${db_dir}"
        publishDir path: { params.saveKrakenDB ? "${params.outdir}/kraken_database" : params.outdir },
            saveAs: { params.saveKrakenDB ? it : null }, mode: 'copy'

        input:
            path db_dir

        output:
            path db_dir

        script:
            download_library = 'download-library'
            build_db = 'build'

            if (params.krn2_task == 'download-taxonomy'){
            """
            build_kraken_db.py --task ${params.krn2_task} --db ${params.krn2_db} --library ${params.krn2_library} -t ${task.cpus} -v
            build_kraken_db.py --task ${download_library} --db ${params.krn2_db} --library ${params.krn2_library} -t ${task.cpus} -v
            build_kraken_db.py --task ${build_db} --db ${params.krn2_db} -t ${task.cpus} -v
            """
            } else if (params.krn2_task == 'standard') {
            """
            build_kraken_db.py --task ${params.krn2_task} --db ${params.krn2_db} -t ${task.cpus} -v
            """
            } 
    }

    process classify_reads {
        label 'high_memory'
        tag "${name}"
        publishDir "${params.outdir}/classification", mode: 'copy',
            saveAs: { filename -> 
            if (filename.indexOf("_report.txt") > 0 ) "report/$filename"
            else if (filename.indexOf("_results.txt") > 0 ) "results/$filename"
            else null
            }

        input:
            tuple val(name), file(reads)
            path db_dir

        output:
            tuple val(name), path("${name}_report.txt"), emit: ch_kraken_report
            tuple val(name), path("${name}_results.txt"), emit: ch_kraken_results

        script:
            if (params.singleEnd){
                if (hasExtension(reads[0], 'gz')) {
                    """
                    kraken2 \\
                        --use-names \\
                        --db ${params.krn2_db} \\
                        --threads ${task.cpus} \\
                        $reads \\
                        --gzip-compressed \\
                        --use-mpa-style \\
                        --report ${name}_report.txt \\
                        --output ${name}_results.txt
                    """
                } else if (hasExtension(reads[1], 'bzip')) {
                    """
                    kraken2 \\
                        --use-names \\
                        --db ${params.krn2_db} \\
                        --threads ${task.cpus} \\
                        $reads \\
                        --bzip2-compressed \\
                        --use-mpa-style \\
                        --report ${name}_report.txt \\
                        --output ${name}_results.txt
                    """
                } else {
                    """
                    kraken2 \\
                        --use-names \\
                        --db ${params.krn2_db} \\
                        --threads ${task.cpus} \\
                        $reads \\
                        --use-mpa-style \\
                        --report ${name}_report.txt \\
                        --output ${name}_results.txt
                    """
                }
            } else {
                if (hasExtension(reads[0], 'gz')) {
                    """
                    kraken2 \\
                        --use-names \\
                        --db ${params.krn2_db} \\
                        --threads ${task.cpus} \\
                        $reads \\
                        --gzip-compressed \\
                        --use-mpa-style \\
                        --report ${name}_report.txt \\
                        --output ${name}_results.txt
                    """
                } else if (hasExtension(reads[0], 'bzip')) {
                    """
                    kraken2 \\
                        --use-names \\
                        --db ${params.krn2_db} \\
                        --threads ${task.cpus} \\
                        $reads \\
                        --bzip2-compressed \\
                        --use-mpa-style \\
                        --report ${name}_report.txt \\
                        --output ${name}_results.txt
                    """
                } else {
                    """
                    kraken2 \\
                        --use-names \\
                        --db ${params.krn2_db} \\
                        --threads ${task.cpus} \\
                        $reads \\
                        --use-mpa-style \\
                        --report ${name}_report.txt \\
                        --output ${name}_results.txt
                    """
                }
            }
    }
    process merge_kraken_reports {
        label 'low_memory'
        tag 'merge_kraken_report'
        publishDir "${params.outdir}/classification", mode: 'copy',
            saveAs: { filename -> 
            if (filename.indexOf(".txt") > 0 ) "merged_tables/$filename"
            else null
            }

        input:
            file reports

        output:
            path "merged_tables.txt", emit: ch_merged_tables
            path "merged_tables_species_summary.txt", emit: ch_merged_species_summary

        script:
            
            """
            merge_mpa_tables.py --input ${reports} -o merged_tables.txt
            """
    }
    process visualize_species_summary {
        label 'low_memory'
        tag 'visualize_species'
        publishDir "${params.outdir}/visualization", mode: 'copy',
            saveAs: {filename -> 
            if (filename.indexOf(".svg") > 0 ) "plots/$filename"
            }
        
        input:
            file (counts)
            file (species_summary)

        output:
            file "*.svg"

        script:
            """
            visualize_abundance.R ${counts} ${species_summary}
            """
    }
}

// check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}