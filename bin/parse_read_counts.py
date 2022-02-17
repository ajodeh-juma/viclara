#!/usr/bin/env python3

"""
parse read counts files
"""

import os
import sys
import argparse
import pandas as pd
from textwrap import dedent


def parse_args():
    """
    """

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(prog="parse_align_stats.py",
                                     formatter_class=CustomFormatter,
                                     argument_default=argparse.SUPPRESS, prefix_chars='--',
                                     description=dedent('''Parse sorted alignment stats''')
                                     )
    required_group = parser.add_argument_group(dedent('''input options'''))
    required_group.add_argument('--stats', metavar="<FILE>", type=str, dest='stats', default=None, help='path to the read count stats file')
    parser.add_argument('--output', metavar='<FILE>', type=str, dest='output', help='output file to write stats')
    return parser


def get_read_stats(stats, output):
    """
    outputs summarized reads stats
    :param stats: path to the stats/log file
    :param output output filename
    """

    stats_dic = dict()
    sampleid = os.path.basename(stats)
    sampleid = os.path.splitext(sampleid)[0]
    print(sampleid)


    for line in open(stats):
        line = line.strip()
        if 'Reads Processed' in line:
            reads = line.strip().split()[2]
            stats_dic[sampleid] = reads

    col_names = {'index': 'Sample', 0: 'Reads'}
    df = pd.DataFrame.from_dict(stats_dic, orient='index').reset_index().rename(columns=col_names)
    print(df)
    df.to_csv(output, index=False, sep="\t")
    return output


def main():
    parser = parse_args()
    args = parser.parse_args()
    get_read_stats(stats=args.stats, output=args.output)

if __name__ == '__main__':
    main()
