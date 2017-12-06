#!/usr/bin/env python3

__author__ = "Sergey Knyazev"
__email__ = "sergey.n.knyazev@gmail.com"
__date__ = "12/6/17"


import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse


AMBIGUOUS = {'AG': 'R',
             'CT': 'Y',
             'GT': 'K',
             'AC': 'M',
             'CG': 'S',
             'AT': 'W',
             'CGT': 'B',
             'AGT': 'D',
             'ACT': 'H',
             'ACG': 'V',
             'ACGT': 'N'
             }

NUCLS = ['A', 'C', 'G', 'T']


def get_nucl(count, thr):
    total = 0.0
    nucls = list()
    for key in NUCLS:
        total += count[key]
    for key in NUCLS:
        if count[key]/total > thr:
            nucls.append(key)
    if len(nucls) >= 2:
        return AMBIGUOUS[''.join(nucls)]
    else:
        return nucls[0]


def get_consensus(file_name, thr):
    samfile = pysam.AlignmentFile(file_name, 'rb')
    consensus = list()
    for col in samfile.pileup():
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        for read in col.pileups:
            if not read.is_del and not read.is_refskip:
                count[read.alignment.query_sequence[read.query_position]] += 1
        consensus.append(get_nucl(count, thr))
    samfile.close()
    return ''.join(consensus)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_bam')
    parser.add_argument('out_fasta')
    parser.add_argument('-t', dest='thr', default=0.1)
    return parser.parse_args()


def main(args):
    SeqIO.write(
        [SeqRecord(Seq(get_consensus(args.in_bam, args.thr)), id='consensus', description='consensus')],
        open(args.out_fasta, 'w'), 'fasta')


if __name__ == '__main__':
    main(parse_args())
