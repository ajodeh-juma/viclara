// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

// check file extension
def checkExtension(file, extension) {
    file.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def checkFile(filePath, extension) {
  // first let's check if has the correct extension
  if (!checkExtension(filePath, extension)) exit 1, "File: ${filePath} has the wrong extension. See --help for more information"
  // then we check if the file exists
  if (!file(filePath).exists()) exit 1, "Missing file in TSV file: ${filePath}, see --help for more information"
  // if none of the above has thrown an error, return the file
  return(file(filePath))
}

// the function expects a comma-separated sample sheet, with a header in the first line
// the header will name the variables and therefore there are a few mandatory names
// id to indicate the sample name or unique identifier
// fastq_1 to indicate fastq_1.fastq.gz, i.e. R1 read or forward read
// fastq_2 to indicate fastq_2.fastq.gz, i.e. R2 read or reverse read
// any other column should fulfill the requirements of modules imported in main
// the function also expects a boolean for single or paired end reads from params

def readInputFile(tsvFile, single_end) {
    Channel.from(tsvFile)
        .splitCsv(header:true, sep: ',')
        .map { row ->
            def meta = [:]
            def reads = []
            def sampleinfo = []
            meta.id = row.sample
            if (single_end) {
              reads = checkFile(row.fastq_1, "fastq.gz")
            } else {
              reads = [ checkFile(row.fastq_1, "fastq.gz"), checkFile(row.fastq_2, "fastq.gz") ]
            }
            sampleinfo = [ meta, reads ]
            return sampleinfo
        }
}