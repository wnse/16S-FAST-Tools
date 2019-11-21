#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 00:21:48 2019

@author: yk
"""

import argparse
import os
import pandas as pd
from Bio import SeqIO

def merge_contigs(contig_list,merge_fa_name):
    seq_no_dict = {}
    with open (merge_fa_name,'w') as w:
        for record in contig_list:
            name = record[0]
            fa = record[1]
            tmp = []
            n = 0
            for seq in SeqIO.parse(fa,'fasta'):
                seq.id = name + '_' + str(n)
                tmp.append(seq)
                n += 1
            seq_no = SeqIO.write(tmp,w,'fasta')
            seq_no_dict[name] = seq_no
    return pd.DataFrame.from_dict(seq_no_dict,orient='index',
                                  columns=(['NO._of_Contigs']))

if __name__ == '__main__':
    parse = argparse.ArgumentParser(description='merge contig fastas in a dir')
    parse.add_argument('-d','--dir',required=True,
                       help='input dir.input dir-dir-fasta')
    parse.add_argument('-o','--outfasta',default='merged.fasta',
                       help='output merge.fasta')
    parse.add_argument('-i','--info',default='merge.fasta.info',
                       help='output number of contigs info of each sample')
    args = parse.parse_args()
    
    input_list = []
    for d in os.listdir(spade_dir):
        contig_file = os.path.join(spade_dir,d,'contigs.fasta')
        if os.path.exists(contig_file):
            tmp = [d,contig_file]
            input_list.append(tmp)
        
    df = merge_contigs(input_list,merge_fa)
    df.to_csv(merge_fa_info,sep='\t')
    