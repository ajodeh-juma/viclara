#!/usr/bin/env python3

'''rename FASTA headers using a specified name (prefix) in masked consensus FASTA file'''

import os
import argparse
import fileinput
from textwrap import dedent


def parse_args():
    """command line options"""
    parser = argparse.ArgumentParser(
        prefix_chars='-',
        description=__doc__
    )

    required_group = parser.add_argument_group(dedent('''mandatory arguments'''))
    required_group.add_argument('-m', '--masked-fasta', metavar='<FILE>', required=True,
                                dest='masked_fasta',
                                help="path to the masked consensus FASTA file"
                                )

    required_group.add_argument('-p', '--prefix', metavar='str', required=True,
                                dest='prefix',
                                help="prefix of the new sample's header names"
                                )

    return parser


def rename_headers(masked_fasta, prefix):
    """
    param masked_fasta
    param prefix
    return
    """

    bak = os.path.join(os.path.dirname(masked_fasta), os.path.basename(masked_fasta) + '.bak')
    if os.path.exists(bak):
        print("renamed masked consensus file exists")
        pass
    else:
        try:
            print('renaming headers/chromosome name(s) in: {}'.format(masked_fasta))
            with fileinput.FileInput(masked_fasta, inplace=True, backup='.bak') as fn:
                for line in fn:
                    if line.startswith('>'):
                        # print('>' + str(prefix) + '_' + line.strip('>'), end='')
                        print('>' + str(prefix), end='\n')
                    else:
                        print(line, end='')
        except Exception as error:
            print(f"Error {error} occurred")
    return masked_fasta


def main():
    parser = parse_args()
    args = parser.parse_args()

    rename_headers(masked_fasta=args.masked_fasta, prefix=args.prefix)


if __name__ == '__main__':
    main()
