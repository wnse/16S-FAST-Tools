#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:26:48 2019

@author: yk
"""

import argparse
import re
from Bio import SeqIO
import pandas as pd

def cut_fa_by_len(input_fa, output_fa, min_len=0, max_len=0):
    name = {}
    name_track = {}
    tmp = []
    if max_len == 0:
        for seq in SeqIO.parse(input_fa, 'fasta'):
            n = re.match('(.+)_\d+', str(seq.id)).group(1)
            l = len(seq.seq)
            if l >= min_len:
                name[n] = name.get(n,0)
                seq.id = n + '_' + str(name[n])
                name_track[seq.id] = name[n]
                name[n] += 1
                tmp.append(seq)
    else:
        for seq in SeqIO.parse(input_fa, 'fasta'):
            n = re.match('(.+)_\d+', seq.id).group(1)
            l = len(seq.seq)
            if l >= min_len and l < max_len:
                name[n] = name.get(n,0)
                seq.id = n + '_' + str(name[n])
                name_track[seq.id] = name[n]
                name[n] += 1
                tmp.append(seq)
    c = SeqIO.write(tmp, output_fa, 'fasta')
    df = pd.DataFrame.from_dict(name_track,orient='index',
                                columns=(['merge_ID']))
    return c, df
            
            
if __name__ == '__main__':
    parse = argparse.ArgumentParser(
            description='get seqs of specific length in fasta',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-i','--infasta',required=True,help='input fasta file')
    parse.add_argument('-o','--outfasta',required=True,help='output fasta file')
    parse.add_argument('-min','--minLength',default=0,type=int,
                       help='min length of seqs')
    parse.add_argument('-max','--maxLength',default=0,type=int,
                       help='max length of seqs. 0 means no max limits')
    args = parse.parse_args()
    
    c, df = cut_fa_by_len(args.infasta, args.outfasta, 
                          args.minLength, args.maxLength)
    print('total seqs:\t{}'.format(c))
    df.to_csv(args.outfasta + '.id.info',sep='\t')
    