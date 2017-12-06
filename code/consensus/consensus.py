#!/usr/bin/env python3

__author__ = "Sergey Knyazev"
__email__ = "sergey.n.knyazev@gmail.com"
__date__ = "12/6/17"


import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


def get_consensus(file_name):
    samfile = pysam.AlignmentFile(file_name, 'rb')
    consensus = list()
    for col in samfile.pileup():
        count = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
        for read in col.pileups:
            if not read.is_del and not read.is_refskip:
                count[read.alignment.query_sequence[read.query_position]] += 1
        consensus.append(max(count, key=count.get))
    samfile.close()
    return ''.join(consensus)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_bam')
    parser.add_argument('out_fasta')
    return parser.parse_args()


def main(args):
    SeqIO.write(
        [SeqRecord(Seq(get_consensus(args.in_bam)), id='consensus', description='consensus')],
        open(args.out_fasta, 'w'), 'fasta')


if __name__ == '__main__':
    main(parse_args())
