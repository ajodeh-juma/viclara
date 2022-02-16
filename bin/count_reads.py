#!/usr/bin/env python3

"""
count reads
"""

import os
import re
import sys
import argparse
import subprocess
import pandas as pd
from textwrap import dedent
from utils import find_executable


def parse_args():
    """
    """

    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(prog="count_reads.py",
                                     formatter_class=CustomFormatter,
                                     argument_default=argparse.SUPPRESS, prefix_chars='--',
                                     description=dedent('''count reads in fastq file''')
                                     )
    required_group = parser.add_argument_group(dedent('''input options'''))
    required_group.add_argument('--fwd', metavar="<FILE>", type=str, dest='fwd', default=None, help='path to the forward FASTQ file')
    parser.add_argument('--rev', metavar="<FILE>", type=str, dest='rev', default=None, help='path to the reverse FASTQ file')
    parser.add_argument('--output', metavar='<FILE>', type=str, dest='output', help='output file to write stats')
    return parser


def get_read_counts(fastq):
    """
    outputs read counts as integer given a fastq file
    :param fastq: path to the FASTQ file
    """

    # locate the executable
    reformat = find_executable(["reformat.sh"])

    # command to count reads
    call = ["{} in={} -Xmx4G".format(reformat, fastq)]
    cmd = "".join(call)

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout_iterator = iter(proc.stdout.readline, b"")

    summary_dict = dict()
    for line in stdout_iterator:
        line = line.decode('utf-8').rstrip()
        m1 = re.match(r"^java.*$", line.strip())
        m2 = re.match(r"^Input:.*$", line.strip())
        if m1:
            read = os.path.basename(m1.group(0).split()[-2]).strip("in=")
            summary_dict[read] = ''
        if m2:
            read_count = int(m2.group(0).split()[1])
            if read in summary_dict:
                summary_dict[read] = read_count
    return summary_dict

def get_pe_rc(fwd, rev, output):
    """
    count PE reads data given the forward ane reverse FASTQ files
    
    :param fwd
    :param rev
    :param output: output stream to which integer read count value is written

    """

    # get sampleids
    samples_dict = dict()
    ext = "fastq fq".split()

    if rev is None:
        f = get_read_counts(fwd)
        fnames = [k for k in f.keys()]
        for fname in fnames:
            if (ext[0] in fname or ext[1] in fname):
                base = os.path.basename(fname)
                sampleid = base.split('.fastq')[0].split('.fq')[0]
                sampleid = sampleid.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                sampleid = sampleid.replace('.R1', '_R1').replace('.R2', '_R2').replace("_R1", "").replace("_R2", "")
                sampleid = re.sub('_S(\d+)_L\d{3}_\d{3}', '', sampleid)
                sampleid = sampleid.replace(".", "_").replace(" ", "_")
                sampleid = re.sub(r'\_.*', '', sampleid)
                # sampleid = re.sub('_L\d{3}_\d{3}', '', sampleid).replace(".trim", "")

                if not sampleid in samples_dict:
                    samples_dict[sampleid] = [f.get(fnames[0])]
        for k, v in samples_dict.items():
            samples_dict[k] = sum(v)
        df = pd.DataFrame(samples_dict.items(), columns=['Sample', 'Reads'])
        #df = pd.DataFrame.from_dict(samples_dict,orient='index').reset_index().rename(columns={'index': 'Sample', 0: 'read1'})
        df.to_csv(output, index=False, sep="\t")
    else:
        f = get_read_counts(fwd)
        r = get_read_counts(rev)
        fnames = [k for k in f.keys()] + [k for k in r.keys()]
        for fname in fnames:
            if (ext[0] in fname or ext[1] in fname):
                base = os.path.basename(fname)
                sampleid = base.split('.fastq')[0].split('.fq')[0]
                sampleid = sampleid.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                sampleid = sampleid.replace('.R1', '_R1').replace('.R2', '_R2').replace("_R1", "").replace("_R2", "")
                sampleid = re.sub('_L\d{3}_\d{3}', '', sampleid).replace("_1.trim", "").replace("_2.trim", "")

                if not sampleid in samples_dict:
                    samples_dict[sampleid] = [f.get(fnames[0])]
                else:
                    samples_dict[sampleid].append(r.get(fnames[1]))

        for k, v in samples_dict.items():
            samples_dict[k] = sum(v)
        df = pd.DataFrame(samples_dict.items(), columns=['Sample', 'Reads'])
        #df = pd.DataFrame.from_dict(samples_dict,orient='index').reset_index().rename(columns={'index': 'Sample', 0: 'read1', 1: 'read2'})
        df.to_csv(output, index=False, sep="\t")
    return output




def main():
    parser = parse_args()
    args = parser.parse_args()
    get_pe_rc(fwd=args.fwd, rev=args.rev, output=args.output)

if __name__ == '__main__':
    main()
