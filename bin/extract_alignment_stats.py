#!/usr/bin/env python3

"""
count reads
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
    required_group.add_argument('--flagstat', metavar="<FILE>", type=str, dest='flagstat', default=None,
                                help='path to the forward sorted alignment flagstat file')
    parser.add_argument('--output', metavar='<FILE>', type=str, dest='output', help='output file to write stats')
    return parser


def get_mapping_stats(flagstat, output):
    """
    outputs alignment stats
    :param flagstat: path to the flagstat file
    :param output output filename
    """

    map_dict = dict()
    sampleid = os.path.basename(flagstat)
    sampleid = sampleid.rsplit(".", 3)[0]

    
    with open(flagstat) as f:
        lines = [line.rstrip() for line in f]
    
    # print(lines[0].split()[0])
    # print(lines[0].split()[2])

    map_dict[sampleid] = [int(lines[4].split()[0]), int(lines[0].split()[0])]

    # for line in open(flagstat):
    #     line = line.strip()

    #     if 'mapped (' in line:
    #         mapped = int(line.split(" ")[0])
    #         if not sampleid in map_dict:
    #             map_dict[sampleid] = [mapped]
    #     elif 'in total (' in line:
    #         total = int(line.split(" ")[0])
    #         # total = int(line.split(" ")[0])
    #         if not total in map_dict:
    #             map_dict[sampleid] += [total]
    #     else:
    #         continue

    print(map_dict)

    col_names = {'index': 'Sample', 0: 'mapped', 1: 'total_paired'}
    df = pd.DataFrame.from_dict(map_dict, orient='index').reset_index().rename(columns=col_names)
    print(df)
    df.to_csv(output, index=False, sep="\t")
    return output


def main():
    parser = parse_args()
    args = parser.parse_args()
    get_mapping_stats(flagstat=args.flagstat, output=args.output)


if __name__ == '__main__':
    main()
