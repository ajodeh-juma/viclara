#!/usr/bin/env python3
"""
build snpeff database
"""

# --- standard imports ---#
import os
import re
import sys
import logging
import argparse
from textwrap import dedent

# --- project specific imports ---#
from utils import mkdir
from utils import copy_file
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

    parser.add_argument('--snpeff-db', metavar='<dir>',
                                dest='snpeff_db', default='snpeff_db',
                                help="path to the SnpEff data directory"
                                )
    parser.add_argument('--config', metavar='<file>',
                                dest='snpeff_config', default='snpeff.config',
                                help="path to the SnpEff config file"
                                )
    required_group.add_argument('--reference', type=str, metavar='<file>',
                                dest='reference', required=True,
                                help='reference file in FASTA format'
                                )
    required_group.add_argument('--gff', metavar='<gff>',
                                dest='gff',
                                help="file with genomic features coordinates in the reference"
                                )
    return parser


def build_snpeff_db(reference, gff, snpeff_config, snpeff_db):
    """
    build SnpEff database for a reference genome
    
    :param: snpeff_config
    :param snpeff_db:
    :param reference:
    :param gff:
    :return:
    """

    # locate the executable
    snpeff = find_executable(['snpEff'])

    snpeff_db = os.path.abspath(snpeff_db)

    # create SnpEff database
    prefix = os.path.join(os.path.abspath(os.path.dirname(reference)), os.path.splitext(os.path.basename(reference))[0])
    index_base = os.path.basename(prefix)
    print(index_base)
    snpeff_data_dir = os.path.join(snpeff_db, 'data')
    snpeff_genes_dir = os.path.join(snpeff_data_dir, index_base)
    mkdir(snpeff_data_dir)
    mkdir(snpeff_genes_dir)

    # copy the files 
    copy_file(src=gff, dest=os.path.join(snpeff_genes_dir, 'genes.gff'))
    copy_file(src=reference, dest=os.path.join(snpeff_genes_dir, 'sequences.fa'))

    # Add a genome to the configuration file
    snpeff_config = os.path.join(snpeff_db, 'snpeff.config')
    with open(snpeff_config, 'w') as f_obj:
        f_obj.write('{}.genome : {}\n'.format(index_base, index_base))

    # check if db exists and build if not
    db_bin = os.path.join(snpeff_genes_dir, 'snpEffectPredictor.bin')
    if os.path.exists(db_bin):
        logging.critical("SnpEff database exist for {}".format(index_base))
    else:
        # build db
        call = ["{} build -config {} -dataDir {} -gff3 -v {}".format(snpeff, snpeff_config, snpeff_data_dir, index_base)]
        cmd = " ".join(call)
        logging.info("building SnpEFF database: {}".format(gff))
        run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return snpeff_config, snpeff_data_dir



def main():
    """

    :return:
    """
    parser = parse_args()
    args = parser.parse_args()

    # build SnpEff database
    build_snpeff_db(
                    reference=args.reference,
                    gff=args.gff,
                    snpeff_db=args.snpeff_db,
                    snpeff_config=args.snpeff_config
                    )


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