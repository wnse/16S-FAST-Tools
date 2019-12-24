#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:06:49 2019

@author: yk
"""

import argparse
import logging
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
def reverse_complement(seq_list):
    '''序列取反向互补
    参数：
        seq_list: 序列列表（tpye=list）
    返回：
        rc_list: 序列列表（type=list）
    '''
    rc_list=[]
    for seq in seq_list:
        rc_list.append(str(Seq(seq,generic_dna).reverse_complement()))
    return rc_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
                                     'fasta sequence reverse complement')
    parser.add_argument('in_fasta',nargs='?',help='input fasta')
    args = parser.parse_args()
    
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
    if args.in_fasta:
        for seq in SeqIO.parse(args.in_fasta,'fasta'):
            print('>' + seq.id + ' reverse_complement')
            print(reverse_complement([str(seq.seq)])[0])
    else:
        seq = input('input sequences:\n')
        print('output reverse complement')
        print(reverse_complement(seq.split()))
        