#!/usr/bin/env python3

"""
merge qc csv files
"""

import os
import re
import sys
import csv
import argparse
import functools
from textwrap import dedent
from itertools import takewhile
import pandas as pd


def parse_args():
    """
    """

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(prog="merge_dataframes.py",
                                     formatter_class=CustomFormatter,
                                     argument_default=argparse.SUPPRESS, prefix_chars='--',
                                     description=dedent('''merge dataframes''')
                                     )
    required_group = parser.add_argument_group(dedent('''input options'''))
    required_group.add_argument('--raw', metavar="<FILE1>", required=True,
                                dest='raw',
                                help='path to the merged raw read counts'
                                )
    required_group.add_argument('--trimmed', metavar="<FILE1>", required=True,
                                dest='trimmed',
                                help='path to the merged trimmed read counts'
                                )
    required_group.add_argument('--mapped', metavar="<FILE1>", required=True,
                                dest='mapped',
                                help='path to the merged mapped read counts'
                                )
    required_group.add_argument('-p', '--prefix', metavar='<STR>', dest='prefix', help='prefix to the output file')
    return parser


def merge_dfs(raw, trimmed, mapped, prefix):
    """
    merge multiple dataframes on common columns

    :param raw:
    :param trimmed:
    :param mapped:
    :return:
    """

    raw = pd.read_csv(raw)
    trimmed = pd.read_csv(trimmed)
    mapped = pd.read_csv(mapped)

    df = functools.reduce(lambda left, right: pd.merge(left, right, on=['Sample'], how='inner'),
                          [raw, trimmed])
    df = functools.reduce(lambda left, right: pd.merge(left, right, on='Sample', how='inner'),
                          [df, mapped])

    df = df.rename(columns={'Reads_x': 'before_trim', 'Reads_y': 'after_trim'})
    # write to output
    output = prefix + '.read.counts.csv'
    df.to_csv(output, index=False, sep=",")
    return output


def main():
    parser = parse_args()
    args = parser.parse_args()
    try:
        merge_dfs(raw=args.raw, trimmed=args.trimmed, mapped=args.mapped, prefix=args.prefix)
    except Exception as e:
        print("I/O Error {}".format(e))


if __name__ == '__main__':
    main()
