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
    required_group.add_argument('--flagstat', metavar="<FILE>", type=str, dest='flagstat', default=None, help='path to the forward sorted alignment flagstat file')
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

    for line in open(flagstat):
        line = line.strip()
        
        if 'mapped (' in line:
            mapped = int(line.split(" ")[0])
            if not sampleid in map_dict:
                map_dict[sampleid] = [mapped]
        elif 'paired in sequencing' in line:
            total = int(line.split(" ")[0])
            if sampleid in map_dict:
                map_dict[sampleid] += [total]
        else:
            continue
            
            
            # perc_mapped = "{0:.2f}".format(perc_mapped)
            # if sampleid in map_dict:
            #     map_dict[sampleid] += [mapped, perc_mapped]
        # elif 'properly paired' in line:
        #     proper_pair = int(line.split(" ")[0])
        #     perc_pp = float(proper_pair/total)*100
        #     perc_pp = "{0:.2f}".format(perc_pp)
        #     if sampleid in map_dict:
        #         map_dict[sampleid] += [proper_pair, perc_pp]

    # col_names = {'index': 'Sample', 0: 'total', 1: 'mapped', 2: 'percentage_mapped', 3: 'proper_pairs', 4: 'percentage_proper_pairs'}

    col_names = {'index': 'Sample', 0: 'mapped', 1: 'total_paired'}
    df = pd.DataFrame.from_dict(map_dict, orient='index').reset_index().rename(columns=col_names)
    df.to_csv(output, index=False, sep="\t")
    return output


def main():
    parser = parse_args()
    args = parser.parse_args()
    get_mapping_stats(flagstat=args.flagstat, output=args.output)

if __name__ == '__main__':
    main()
