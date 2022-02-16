#!/usr/bin/env python3

"""
merge qc csv files
"""

import os
import re
import sys
import csv
import argparse
from textwrap import dedent
from itertools import takewhile
import pandas as pd


def parse_args():
    """
    """

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(prog="merge_read_counts.py",
                                     formatter_class=CustomFormatter,
                                     argument_default=argparse.SUPPRESS, prefix_chars='--',
                                     description=dedent('''merge read counts''')
                                     )
    required_group = parser.add_argument_group(dedent('''input options'''))
    required_group.add_argument('--input', metavar="<FILE1,FILE2.....FILEn>", nargs='+', required=True,
                                dest='input',
                                help='one or more tab-separated files to join, paths whould be separated by space'
                                )
    required_group.add_argument('-p', '--prefix', metavar='<STR>', dest='prefix', help='prefix to the output file')
    required_group.add_argument('-f', '--ftype', metavar="<STR>", dest="ftype", help="count category", 
                                choices=['raw', 'trimmed', 'mapped'])
    return parser


def merge_dfs(input, prefix, ftype):
    """
    outputs the table join of the given pre-split string collection
    :param input: collection of collections of string collections
    :param prefix: prefix to the output file
    :param ftype file type as either read-count or flagstat
    """

    csv_dict = dict()    
    # reads_bool = [f.endswith('.tsv') for f in input]
    # flagstat_bool = [f.endswith('.flagstat') for f in input]
    # rc_res, fl_res = all(reads_bool), all(flagstat_bool)
    
    
    for f in input:
        sample = os.path.splitext(os.path.basename(f))[0].strip('.tsv')
        with open(f) as fin:
            for line in fin:
                if line.startswith('Sample'):
                    header = line.strip().split('\t')
                    continue
                else:
                    csv_dict[sample] = line.strip().split('\t')
    df = pd.DataFrame.from_dict(csv_dict, orient='index', columns=header)

    # write to output
    if ftype == 'raw':
        output = prefix + '.raw.csv'
    elif ftype == 'trimmed':
        output = prefix + '.trimmed.csv'
    else:
        output = prefix + '.mapped.csv'
    df.to_csv(output, index=False, sep=",")


def main():
    parser = parse_args()
    args = parser.parse_args()
    try:
        merge_dfs(input=args.input, prefix=args.prefix, ftype=args.ftype)
    except Exception as e:
        print("I/O Error {}".format(e))


if __name__ == '__main__':
    main()
