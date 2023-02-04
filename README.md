<!-- # ![nf-core/viclara](docs/images/nf-core-viclara_logo.png) -->

<!-- **Identify Rift Valley Fever virus lineages of a nucleotide sequence**. -->

  - [Introduction](#introduction)
  - [Installation](#installation)
  - [Testing](#test)
  - [Usage](#usage)
  - [Database](#database)
  - [Method details](#method-details)
  - [Output](#output)
  - [Pipeline summary](#pipeline-summary)
  - [Citations](#citations)
  

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23viclara-4A154B?logo=slack)](https://nfcore.slack.com/channels/viclara)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/ajodeh-juma/viclara/blob/master/LICENSE)
<!-- [![Docker](https://img.shields.io/docker/automated/nfcore/viclara.svg)](https://hub.docker.com/r/nfcore/viclara) -->
<!-- [![GitHub Actions CI Status](https://github.com/nf-core/viclara/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/viclara/actions) -->
<!-- [![GitHub Actions Linting Status](https://github.com/nf-core/viclara/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/viclara/actions) -->
<!-- ![](https://img.shields.io/badge/uses-docker-blue.svg) -->

[![Twitter Follow](https://img.shields.io/twitter/follow/john_juma.svg?style=social)](https://twitter.com/john_juma)

## Introduction

**ViCLaRA** is a bioinformatics analysis pipeline for classification and reference guided assembly of segmented viruses from metagenomics reads obtained on Illumina platform. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.


## Installation

**ViCLaRA** runs on UNIX/LINUX systems. You will install Miniconda3 from [here](https://docs.conda.io/en/latest/miniconda.html). Once Miniconda3 has been installed, proceed with pipeline installation

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline

```
git clone https://github.com/ajodeh-juma/viclara.git
cd rvfvMETA
conda env create -n viclara-dev1.0 -f environment.yml
conda activate viclara-dev1.0
```

## Testing

  - Optional: Test the installation on a minimal dataset bundled in the installation
    
    ```nextflow run main.nf -profile test```
    

## Usage

For minimal pipeline options, use the ```--help``` flag e.g. 

```nextflow run main.nf --help```

To see all the options, use the ```--show_hidden_params``` flag e.g.

```nextflow run main.nf --help --show_hidden_params```

A typical command for S segment is
```
nextflow run main.nf \
   --input 'data/test/*_R{1,2}.fastq.gz' \
   --segment S \
   --outdir output-dir \
   -work-dir work-dir \
```

## Database

If you intend to identify off-target reads in your samples, the pipeline uses `Kraken2` to classify the reads.
One can create a standard Kraken2 database (which is quite large and requires ample memory when creating ) or a
custom database with organisms of interest. Once created, the database can be given as an option while running the pipeline

Typically you can create a database using the command:
```
# install the taxonony
DBNAME="/Users/jjuma/databases/kraken2/viral/"
kraken2-build --download-taxonomy --db $DBNAME

# install a reference library. In this case, `viral` which has all RefSeq complete viral genomes/proteins
kraken2-build --download-library viral --db $DBNAME

# build the database
kraken2-build --build --db $DBNAME

```

You can then run the pipeline as
```
nextflow run main.nf \
   --input 'data/test/*_R{1,2}.fastq.gz' \
   --segment S \
   --database /Users/jjuma/databases/kraken2/viral/ \
   --outdir output-dir \
   -work-dir work-dir \
```

## Method details

The pipeline offers several parameters including as highlighted:

```
Input/output options
  --input                      [string]  Input FastQ files.
  --single_end                 [boolean] Specifies that the input is single-end reads.
  --metadata                   [string]  Input comma-separated values (csv) metadata file containing the columns 'sample_name' and  'Ct.
  --outdir                     [string]  The output directory where the results will be saved. [default: ./results]
  --email                      [string]  Email address for completion summary.

genome options
  --segment                    [string]  genomic segment of the virus. options are 'S', 'M' and 'L'
  --filter_phix                [boolean] Filter PhiX sequences
  --host_fasta                 [string]  Path to the host FASTA genome file. Not required if either --host_bwa_index or --host_bowtie2_index is 
                                         specified 
  --host_bwa_index             [string]  Path to host genome directory or tar.gz archive for pre-built BWA index.
  --host_bowtie2_index         [string]  Path to host genome directory or tar.gz archive for pre-built BOWTIE2 index.

Read trimming options
  --trimmer                    [string]  Specifies the alignment algorithm to use - available options are 'fastp', 'trimmomatic'. [default: fastp]
  --adapters                   [string]  Path to FASTA adapters file
  --illumina_clip              [string]  Instructs Trimmomatic apply a maximum mismatch count, how accurate the match between the two adapter ligated reads 
                                         reads must be for PE palindrome, and how accurate the match between any adapter etc. sequence must be against a 
                                         read [default: 2:30:10] 
  --leading                    [integer] Instructs Trimmomatic to cut bases off the start of a read, if below a threshold quality [default: 3]
  --trailing                   [integer] Instructs Trimmomatic to cut bases off the end of a read, if below a threshold quality [default: 3]
  --average_quality            [integer] Instructs Trimmomatic or Fastp the average quality required in the sliding window [default: 20]
  --min_length                 [integer] Instructs Trimmomatic or Fastp to drop the read if it is below a specified length [default: 50]
  --window_size                [integer] Instructs Trimmomatic to apply a sliding window (number of bases to average across) [default: 4]
  --window_quality             [integer] Instructs Trimmomatic to apply the average quality required in the sliding window [default: 20]
  --qualified_quality_phred    [integer] Instructs Fastp to apply the --qualified_quality_phred option [default: 30]
  --unqualified_percent_limit  [integer] Instructs Fastp to apply the --unqualified_percent_limit option [default: 10]
  --skip_trimming              [boolean] Skip the adapter trimming step.
  --save_trimmed               [boolean] Save the trimmed FastQ files in the results directory.
  --save_trimmed_fail          [boolean] Save failed trimmed reads.

Alignment options
  --aligner                    [string]  Specifies the alignment algorithm to use - available options are 'bwa', 'bowtie2'. [default: bwa]
  --seq_center                 [string]  Sequencing center information to be added to read group of BAM files.
  --save_align_intermeds       [boolean] Save the intermediate BAM files from the alignment step.
  --min_mapped                 [integer] Minimum number of mapped reads to be used as threshold to drop low mapped samples [default: 200]

Variant calling options
  --variant_caller             [string]  Specifies the variant calling algorithm to use - available options are 'bcftools', 'varscan'. [default: 
                                         bcftools] 
  --mpileup_depth              [integer] SAMtools mpileup max per-file depth, avoids excessive memory usage
  --min_base_quality           [integer] Skip bases with baseQ/BAQ smaller than this value when performing variant calling
  --min_coverage               [integer] Skip positions with an overall read depth smaller than this value when performing variant calling
  --min_allele_freq            [number]  Minimum allele frequency threshold for calling variant  [default: 0.25]
  --max_allele_freq            [number]  Maximum allele frequency threshold for calling variant  [default: 0.75]
  --varscan_strand_filter      [boolean] Ignore variants with >90% support on one strand [default: true]
  --save_mplieup               [boolean] Save SAMtools mpileup output file

Process skipping options
  --skip_multiqc               [boolean] Skip MultiQC.
  --skip_qc                    [boolean] Skip all QC steps except for MultiQC.
  --skip_trimming              [boolean] Skip quality trimming step.
  --skip_kraken2               [boolean] Skip classification of reads using Kraken2.
  --skip_markduplicates        [boolean] Skip picard MarkDuplicates step.
  --skip_alignment             [boolean] Skip all of the alignment-based processes within the pipeline.
  --skip_assembly              [boolean] Skip De novo assembly step.

Reads classification options
  --database                   [string]  Path to the kraken2 database name
```

## Output

All the output results will be written to the `results` directory if no `--outdir` is not used. Masked and non-masked consensus genomes will be located in `bcftools/consensus`


See [usage docs](https://github.com/ajodeh-juma/viclara/usage) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:

* Sequencing quality control (`FastQC`)
* Reads classification (`Kraken2`)
* Quality control and preprocessing (`fastp`) or (`trimmomatic`)
* Reads alignment/mapping (`BWA`) or (`Bowtie2`)
* Alignment summary (`SAMtools`)
* Call variants (`BCFtools`) or (`Varscan2`)
* Annotate variants (`SnpEff`) or (`SnpSift`)
* Genome coverage (`BEDTools`)
* Visualization (`R`), (`ggplot2`)
* De novo assembly (`MetaSPAdes`) or (`Trinity`) or (`Tadpole`)
* Overall pipeline run summaries (`MultiQC`)

## Documentation
Generating whole genome sequences of segmented viruses has largely depended on sequencing of partial gene sequences of the viruses. Here we implement a pipeline
that can be adopted to other segmented viruses in order to assemble complete genomic sequences from RNA metagenomic sequencing. We implement this pipeline to
generate complete genome sequences of Rift Valley fever virus, a tripartite virus having 3 segments - Small (S), Medium (M) and Large (L).

The pipeline comes bundled with reference genome and annotation, and the user only has to specify the segment to obtain full genome sequences. The pipeline calls variants
using 2 commonly used tools: `bcftools` or `varscan` and annotates the variants using `SnpEff` and `SnpSift`

## Credits

**ViCLaRA** was originally written by @ajodeh-juma with inspiration from the @nf-core team, particularly on viralrecon

We thank the following people for their extensive assistance in the development
of this pipeline:

## License
ViCLaRA is free software, licensed under [GPLv3](https://github.com/ajodeh-juma/viclara/blob/master/LICENSE).

## Issues
Please report any issues to the [issues page](https://github.com/ajodeh-juma/viclara/issues).

## Contribute
If you wish to fix a bug or add new features to the software we welcome Pull Requests. We use
[GitHub Flow style development](https://guides.github.com/introduction/flow/). Please fork the repo, make the change, then submit a Pull Request against out master branch, with details about what the change is and what it fixes/adds. 
We will then review your changes and merge them, or provide feedback on enhancements.


## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/viclara for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

In addition, references of tools and data used in this pipeline are as follows:
> Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data. http://www.bioinformatics.babraham.ac.uk/projects/fastqc
>
> Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. 
> _Bioinformatics_, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560
>
> Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. 
> _Bioinformatics (Oxford, England)_, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170
>
> Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. 
> _Bioinformatics_, 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324
>
> Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. _Nat Methods_. 2012;9(4):357-359. Published 2012 Mar 4. doi:10.1038/nmeth.1923
>
> Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). 
> The Sequence Alignment/Map format and SAMtools. 
> _Bioinformatics_, 25(16)> 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
>
> Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. 
> _Bioinformatics_. 2011;27(21):2987-2993. doi:10.1093/bioinformatics/btr509
>
> Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, Wilson RK. 
> VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. 
> _Genome Res_. 2012 Mar;22(3):568-76. doi: 10.1101/gr.129684.111. Epub 2012 Feb 2. PMID: 22300766; PMCID: PMC3290792.
>
> R Core Team. (2017). R: A language and environment for statistical computing. _R Foundation for Statistical Computing_,. https://www.R-project.org/
>
> Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. _Genome Biology_, 20(1), 257. https://doi.org/10.1186/s13059-019-1891-0
>
> Cingolani, P., Platts, A., Wang, l., Coon, M., Nguyen, T., Wang, L., Land, S. J., Lu, X., & Ruden, D. M. (2012). 
> A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. 
> _Fly_, 6(2), 80–92. https://doi.org/10.4161/fly.19695