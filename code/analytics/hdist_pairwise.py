#!/usr/bin/env python3

__author__ = "Sergey Knyazev"
__email__ = "sergey.n.knyazev@gmail.com"
__date__ = "12/6/17"


import argparse
from Bio import SeqIO

AMBIGUOUS = {'R': ['A', 'G'],
             'Y': ['C', 'T'],
             'K': ['G', 'T'],
             'M': ['A', 'C'],
             'S': ['C', 'G'],
             'W': ['A', 'T'],
             'B': ['C', 'G', 'T'],
             'D': ['A', 'G', 'T'],
             'H': ['A', 'C', 'T'],
             'V': ['A', 'C', 'G'],
             'N': ['A', 'C', 'G', 'T']
             }


def nucl_dist(a, b):
    if a in AMBIGUOUS:
        a = set(AMBIGUOUS[a])
    else:
        a = {a}
    if b in AMBIGUOUS:
        b = set(AMBIGUOUS[b])
    else:
        b = {b}
    return not a.intersection(b)


def h_dist(seq1, seq2):
    return sum(map(lambda x: nucl_dist(x[0], x[1]), zip(seq1, seq2)))


def print_h_dist(fastas):
    seqs = list()
    for f in fastas:
        seqs += list(SeqIO.parse(open(f), 'fasta'))
    for i in range(len(seqs)):
        print(seqs[i].id, end='\t')
        for j in range(len(seqs)):
            print(h_dist(seqs[i].seq, seqs[j].seq), end='\t')
        print()


def parse_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('fasta', nargs='+')
    return arg_parser.parse_args()


def main(args):
    print_h_dist(args.fasta)


if __name__ == '__main__':
    main(parse_args())
