#!/usr/bin/env python3
"""
Annotate variants using SnpEff

"""

# --- standard imports ---#
import os
import re
import sys
import logging
import argparse
from textwrap import dedent

# --- third party imports ---#


# --- project specific imports ---#
from utils import mkdir
from utils import copy_file
#from utils import nthreads
#from utils import ref_strings
from utils import find_executable
from utils import run_shell_command

log_level = logging.DEBUG
logging.basicConfig(level=log_level,
                    format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p'
                    )


def parse_args():
    """command line options"""

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        formatter_class=CustomFormatter,
        prefix_chars='-',
    )

    required_group = parser.add_argument_group(dedent('''required arguments'''))

    required_group.add_argument('-v', '--vcf-file', type=str, metavar='<dir>',
                                dest='vcf_file',
                                help="path to the variants file (vcf file format)"
                                )

    required_group.add_argument('-db', '--snpeff-db', metavar='<dir>',
                                dest='snpeff_db', default='snpeff_db',
                                help="path to the SnpEff data directory"
                                )
    parser.add_argument('-c', '--config', metavar='<file>',
                        dest='snpeff_config', default='snpeff.config',
                        help="path to the SnpEff config file"
                        )
    required_group.add_argument('-r', '--reference', type=str, metavar='<file>',
                                dest='reference', required=True,
                                help='reference file in FASTA format'
                                )
    required_group.add_argument('-p', '--prefix', type=str, metavar='<str>',
                                dest='prefix', required=True,
                                help='prefix to the output file'
                                )
    required_group.add_argument('-g', '--gff', metavar='<gff>',
                                dest='gff',
                                help="file with genomic features coordinates in the reference"
                                )
    return parser


def build_snpeff_db(reference, gff, snpeff_db, snpeff_config):
    """
    build SnpEff database for a reference genome
    
    :param reference:
    :param gff:
    :param snpeff_db:
    :return:
    """

    # locate the executable
    snpeff = find_executable(['snpEff'])

    snpeff_db = os.path.abspath(snpeff_db)

    # create SnpEff database
    prefix = os.path.join(os.path.abspath(os.path.dirname(reference)), os.path.splitext(os.path.basename(reference))[0])
    index_base = os.path.basename(prefix)
    snpeff_data_dir = os.path.join(snpeff_db, 'data')
    snpeff_genes_dir = os.path.join(snpeff_data_dir, index_base)
    mkdir(snpeff_data_dir)
    mkdir(snpeff_genes_dir)

    # copy the files 
    copy_file(src=gff, dest=os.path.join(snpeff_genes_dir, 'genes.gff'))
    # copy_file(src=gff, dest=os.path.join(snpeff_genes_dir, 'genes.gtf'))
    copy_file(src=reference, dest=os.path.join(snpeff_genes_dir, 'sequences.fa'))

    # Add a genome to the configuration file
    snpeff_config = os.path.join(snpeff_db, snpeff_config)
    with open(snpeff_config, 'w') as f_obj:
        f_obj.write('{}.genome : {}\n'.format(index_base, index_base))

    # check if db exists and build if not
    db_bin = os.path.join(snpeff_genes_dir, 'snpEffectPredictor.bin')
    if os.path.exists(db_bin):
        logging.critical("SnpEff database exist for {}".format(index_base))
    else:
        # build db
        call = ["{} build -config {} -dataDir {} -gff3 -v {}".format(snpeff, snpeff_config, snpeff_data_dir, index_base)]
        # call = ["{} build -config {} -dataDir {} -gtf22 -v {}".format(snpeff, snpeff_config, snpeff_data_dir, index_base)]
        cmd = " ".join(call)
        logging.info("building SnpEFF database: {}".format(gff))
        run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return index_base, snpeff_config, snpeff_data_dir


def annotate_snps(index_base, config, vcf_file, db, snpeff_csv, snpeff_vcf):
    """
    annotate and predict the effect of variants

    :param index_base:
    :param config:
    :param vcf_file:
    :param db:
    :param snpeff_csv:
    :param snpeff_vcf:
    :param out_dir:
    :return:
    """

    # locate the executable
    snpeff = find_executable(['snpEff'])

    if os.path.exists(snpeff_vcf):
        logging.critical("variant annotation file {} exists".format(snpeff_vcf))
    else:
        call = ["snpEff {} -config {} -dataDir {} {} -csvStats {} | bgzip -c > {}"
                "".format(index_base, config, db, vcf_file, snpeff_csv, snpeff_vcf)]
        cmd = " ".join(call)
        logging.info("annotating variants and predicting effects: {}".format(vcf_file))
        p = run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
        if p:
            # index the vcf file
            index_vcf(vcf_file=snpeff_vcf)

            sample = os.path.basename(snpeff_vcf).rsplit(".", 2)[0]
            out_dir = os.path.dirname(snpeff_vcf)

            html = os.path.join(os.getcwd(), 'snpEff_summary.html')
            snpeff_html = os.path.join(out_dir, sample + ".snpEff.summary.html")
            copy_file(src=html, dest=snpeff_html)
            os.remove(html)
    return snpeff_vcf


def index_vcf(vcf_file):
    """
    index VCF file with tabix
    :param vcf_file: <str> VCF file to be indexed
    :return:
    """

    call = ["tabix -p vcf -f {}".format(vcf_file)]
    cmd = " ".join(call)
    logging.info("indexing VCF file")
    run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return vcf_file


def get_stats(vcf_file, stats_file):
    """
    calculate variant call stats with bcftools
    :param vcf_file: <str> VCF file
    :param stats_file: <str> file to write stats
    :return:
    """

    # locate the executable
    bcftools = find_executable(['bcftools'])

    call = ["{} stats {} > {}".format(bcftools, vcf_file, stats_file)]
    cmd = " ".join(call)
    logging.info("calculating stats ")
    run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return stats_file


def filter_variants(vcf_file):
    """

    :param vcf_file:
    :return:
    """

    # locate the executable
    snpsift = find_executable(['SnpSift'])

    vcf_pl = 'vcfEffOnePerLine.pl'
    pwd = os.path.dirname(__file__)
    vcf_path = os.path.join(pwd, vcf_pl)
    print(vcf_path)

    sample = os.path.basename(vcf_file).rsplit(".", 2)[0]
    snpsift_file = os.path.join(os.path.dirname(vcf_file), sample + '.snpSift.table.txt')

    if os.path.exists(snpsift_file):
        logging.critical("SnpSift file {} exists!".format(snpsift_file))
    else:
        call = ['zcat {} | {} | {} extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].IMPACT" '
                '"ANN[*].EFFECT" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" '
                '"ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" '
                '"ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" '
                '"EFF[*].AA_LEN" > {}'.format(vcf_file, vcf_path, snpsift, snpsift_file)]
        cmd = " ".join(call)
        print(cmd)
        # filter
        run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return snpsift_file


def main():
    """

    :return:
    """
    parser = parse_args()
    args, remaining = parser.parse_known_args()
    parser.add_argument('--csv', metavar='<csv file>',
                        dest='snpeff_csv',
                        default=os.path.basename(args.vcf_file).rsplit(".", 2)[0] + '.snpEff.csv',
                        help="comma-separated values format file to write stats"
                        )
    parser.add_argument('--vcf-out', metavar='<vcf file>',
                        dest='snpeff_vcf',
                        default=os.path.basename(args.vcf_file).rsplit(".", 2)[0] + '.snpEff.vcf.gz',
                        help="file to write the output in vcf format (vcf.gz)"
                        )
    args = parser.parse_args()

    print(args)

    # build SnpEff database
    index_base, snpeff_config, snpeff_data_dir = build_snpeff_db(
        reference=args.reference,
        gff=args.gff,
        snpeff_db=args.snpeff_db,
        snpeff_config=args.snpeff_config
    )

    # match = re.match(r'^(^.*\.vcf.gz$)', os.path.basename(args.vcf_file))
    # sample = match.group(0).rsplit(".", 2)[0]
    sample = args.prefix
    vcf_file = args.vcf_file

    snpeff_csv = sample + '.snpEff.csv'
    snpeff_vcf = sample + '.snpEff.vcf.gz'
    snpeff_out = annotate_snps(index_base=index_base,
                               config=snpeff_config,
                               vcf_file=vcf_file,
                               db=snpeff_data_dir,
                               snpeff_csv=snpeff_csv,
                               snpeff_vcf=snpeff_vcf
                               )
    filter_variants(vcf_file=snpeff_out)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        logging.debug('', exc_info=True)
        try:
            sys.exit(e.errno)
        except AttributeError:
            sys.exit(1)
