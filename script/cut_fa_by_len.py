#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:26:48 2019
筛选fasta文件符合长度要求的序列
@author: yk
"""

import argparse
from Bio import SeqIO


def cut_fa_by_len(input_fa, output_fa, min_len=0, max_len=0):
    '''筛选fasta文件符合长度要求的序列'''
    tmp = []
    if max_len == 0:
        for seq in SeqIO.parse(input_fa, 'fasta'):
            l = len(seq.seq)
            if l >= min_len:
                seq.description = 'len=' + str(l)
                tmp.append(seq)
    else:
        for seq in SeqIO.parse(input_fa, 'fasta'):
            l = len(seq.seq)
            if l >= min_len and l < max_len:
                seq.description = 'len=' + str(l)
                tmp.append(seq)
    c = SeqIO.write(tmp, output_fa, 'fasta')
    return c


if __name__ == '__main__':
    parse = argparse.ArgumentParser(
        description='get seqs of specific length in fasta',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-i', '--infasta', required=True, help='input fasta file')
    parse.add_argument('-o', '--outfasta', required=True, help='output fasta file')
    parse.add_argument('-min', '--minLength', default=0, type=int,
                       help='min length of seqs')
    parse.add_argument('-max', '--maxLength', default=0, type=int,
                       help='max length of seqs. 0 means no max limits')
    args = parse.parse_args()
    c = cut_fa_by_len(args.infasta, args.outfasta,
                      args.minLength, args.maxLength)
    print('total seqs:\t{}'.format(c))
