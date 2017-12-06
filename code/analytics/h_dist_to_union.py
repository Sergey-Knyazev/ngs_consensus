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


def dist_to_union(ref, seqs):
    nucl_dists = 0
    for i in range(len(ref.seq)):
        nucl_dist_bool = True
        for j in range(len(seqs)):
            if not nucl_dist(ref[i], seqs[j][i]):
                nucl_dist_bool = False
                continue
        if nucl_dist_bool:
            nucl_dists += 1
    return nucl_dists


def get_dist(seq, fasta):
    seqs = list(SeqIO.parse(open(fasta), 'fasta'))
    seq_ind = None
    for i in range(len(seqs)):
        if seqs[i].id == seq:
            seq_ind = i
            break
    ref = seqs[seq_ind]
    seqs = seqs[:seq_ind] + seqs[seq_ind + 1:]
    return dist_to_union(ref, seqs)


def parse_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('seq_of_interest')
    arg_parser.add_argument('fasta')
    return arg_parser.parse_args()


def main(args):
    print(get_dist(args.seq_of_interest, args.fasta))


if __name__ == '__main__':
    main(parse_args())
