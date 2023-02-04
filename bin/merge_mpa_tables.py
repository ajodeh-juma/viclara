#!/usr/bin/env python3
"""
merge mpa-style kraken tables into a single table with clades and relative abundance
"""

# --- standard python imports ---#

import os
import logging
import argparse
import textwrap
from itertools import takewhile

# --- third-party imports ---#
import pandas as pd

log_level = logging.DEBUG
logging.basicConfig(level=log_level,
                    format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p'
                    )


def parse_args():
    """
    command line options
    :return:
    """

    parser = argparse.ArgumentParser(prog="merge_mpa_tables.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=textwrap.dedent('''\
                                            performs a table join on one or more kraken2 mpa-styled report files
                                            ----------------------------------------------------------------------------

                                            input: list of kraken2 report files in mpa format - files/paths separated by space
                                            prefix: prefix of the output files
                                            '''),
                                     argument_default=argparse.SUPPRESS)
    required_group = parser.add_argument_group('''required arguments''')
    required_group.add_argument("-i", "--input", metavar="<input.txt>", nargs="+", dest="input",
                                help="One or more tab-delimited kraken2 report tables (in mpa format) to join")
    required_group.add_argument('--prefix', type=str, metavar="<str>", dest="prefix",
                                help="prefix of the output files")
    return parser


def merge(input, prefix):
    """
    outputs the table join of the given pre-split string collection.
    :param input: <str> collection of collections of string collections
    :param prefix: <str> output prefix
    :return: output files
    """

    merged_tables = pd.DataFrame()
    names = ['clade_name']

    for f in input:
        if os.path.exists(f) and os.path.getsize(f) == 0:
            continue
        else:
            # sample = os.path.splitext(os.path.basename(f))[0].strip('kraken2.report')
            sample = os.path.splitext(os.path.basename(f))[0].rsplit('.')[0]
            names.append(sample)
            with open(f) as fin:
                headers = [x.strip() for x in takewhile(lambda x: x.startswith('#'), fin)]
                if len(headers) == 0:
                    headers = []
                    index_col = [0]

            df = pd.read_csv(f,
                             sep='\t',
                             skiprows=len(headers),
                             names=names,
                             index_col=index_col
                             )

            # convert the number of reads a percentage
            val = df[names[1]].iloc[0]
            df = df.div(val, level=0) * 100

            if merged_tables.empty:
                merged_tables = df.iloc[:, 0].rename(sample).to_frame()
            else:
                merged_tables = pd.merge(
                    df.iloc[:, 0].rename(sample).to_frame(),
                    merged_tables,
                    how='outer',
                    left_index=True,
                    right_index=True
                )

    output = prefix + '.kraken2.tables.txt'
    merged_tables.fillna('0').reset_index().to_csv(output, index=False, sep='\t')

    # get the species level table
    df = merged_tables.fillna(0).reset_index()
    species_df = df[df.clade_name.str.split('|').str[-1].str.startswith('s__')]
    sp_df = species_df.copy()
    sp_df['clade_name'] = sp_df['clade_name'].str.split('|').str[-1].str.replace('s__', '')
    sp_df = sp_df.rename(columns={'clade_name': 'sample'})

    species_summary = prefix + '.kraken2.species.txt'
    print(species_summary)
    sp_df.to_csv(species_summary, index=False, sep="\t")

    return output


def main():
    parser = parse_args()
    args = parser.parse_args()
    merge(input=args.input, prefix=args.prefix)


if __name__ == '__main__':
    main()