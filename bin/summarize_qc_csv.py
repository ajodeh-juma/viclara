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

    parser = argparse.ArgumentParser(prog="summarize_qc_csv.py",
                                     formatter_class=CustomFormatter,
                                     argument_default=argparse.SUPPRESS, prefix_chars='--',
                                     description=dedent('''merge qc summary files''')
                                     )
    required_group = parser.add_argument_group(dedent('''input options'''))
    required_group.add_argument('--input', metavar="<input.csv>", nargs='+', required=True,
                                dest='input',
                                help='one or more comma-separated files to join, paths whould be separated by space'
                                )
    required_group.add_argument('--read-counts', metavar="<FILE>", type=str, required=True,
                                dest='read_counts',
                                help='File having merged read counts'
                                )
    required_group.add_argument('-p', '--prefix', metavar='<STR>', dest='prefix', help='prefix to the output file')
    return parser


def merge_csv(input, read_counts, prefix):
    """
    outputs the table join of the given pre-split string collection
    :param input: collection of collections of string collections
    :param prefix: output file prefix
    """

    csv_dict = dict()
    for f in input:
        sample = os.path.splitext(os.path.basename(f))[0].strip('.qc.csv')
        with open(f) as fin:
            for line in fin:
                if line.startswith('sample_name'):
                    header = line.strip().split(',')
                    continue
                else:
                    csv_dict[sample] = line.strip().split(',')
    df = pd.DataFrame.from_dict(csv_dict, orient='index', columns=header)

    rc = pd.read_csv(read_counts)
    rc = rc.rename(columns={"Sample": "sample_name"})
    rc['pct_mapped'] = rc['mapped']/rc['after_trim'] * 100
    rc['pct_mapped'] = rc['pct_mapped'].astype(float).round(2)
    merged = pd.merge(rc, df, on='sample_name')
    merged = merged.drop(columns=['total_paired', 'num_aligned_reads'])
    output = prefix + '.quality.control.csv'
    merged.to_csv(output, index=False, sep=",")


def main():
    parser = parse_args()
    args = parser.parse_args()
    try:
        merge_csv(input=args.input, read_counts=args.read_counts, prefix=args.prefix)
    except Exception as e:
        print("I/O Error {}".format(e))


if __name__ == '__main__':
    main()
