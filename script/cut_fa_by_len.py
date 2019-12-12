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
#             n = re.match('(.+)_\d+', str(seq.id)).group(1)
            l = len(seq.seq)
            if l >= min_len:
#                 old_id = seq.id
#                 name[n] = name.get(n,0)
#                 new_id = str(len(name.keys())-1) + '_' + str(name[n])
#                 name_track[old_id] = new_id
#                 seq.id = new_id
#                 name[n] += 1
                seq.description = 'len='+str(l)
                tmp.append(seq)
    else:
        for seq in SeqIO.parse(input_fa, 'fasta'):
#             n = re.match('(.+)_\d+', seq.id).group(1)
            l = len(seq.seq)
            if l >= min_len and l < max_len:
#                 old_id = seq.id
#                 name[n] = name.get(n,0)
#                 new_id = str(len(name.keys())-1) + '_' + str(name[n])
#                 name_track[old_id] = new_id
#                 seq.id = new_id
#                 name[n] += 1
                seq.description = 'len='+str(l)
                tmp.append(seq)
    c = SeqIO.write(tmp, output_fa, 'fasta')
#     df = pd.DataFrame.from_dict(name_track,orient='index',
#                                 columns=(['CutLen_ID']))
#     return c, df
    return c
            
            
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
    
    c = cut_fa_by_len(args.infasta, args.outfasta, 
                          args.minLength, args.maxLength)
    print('total seqs:\t{}'.format(c))
#     df.to_csv(args.outfasta + '.id.info',sep='\t')
    