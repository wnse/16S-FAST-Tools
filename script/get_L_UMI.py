#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:32:45 2019

@author: yk
"""
import argparse
import logging
import pandas as pd
from Bio import SeqIO

def get_L_UMI(LR1,LR2):
    umiL = {}
    for idx in LR1.keys():
        seq1=str(LR1[idx].seq)
        seq2=str(LR2[idx].seq)
        u = ''
        if (len(seq1) != 14):
            seq1 = '*'
        if (len(seq2) != 14):
            seq2 = '*'
        u = seq1 + '|' + seq2
        umiL[u] = umiL.get(u,0)+1
    df = pd.DataFrame.from_dict(umiL,orient='index',columns=(['counts_of_paire'])).reset_index()
    df['u1'],df['u2'] = df['index'].str.split('|').str
    dfu1 = df.groupby(['u1']).sum().rename(columns={'counts_of_paire':'u1_reads'})
    dfu2 = df.groupby(['u2']).sum().rename(columns={'counts_of_paire':'u2_reads'})
    df = pd.merge(df,dfu1,left_on='u1',right_index=True,how='left')
    df = pd.merge(df,dfu2,left_on='u2',right_index=True,how='left')
    #df = df[~df['index'].str.contains('\*')]
    df = df.sort_values(by='counts_of_paire',ascending=False).drop('index',axis=1).reset_index().drop('index',axis=1)
    return df

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
    parse = argparse.ArgumentParser(description='get L data UMI from paired *cut_primer.fq,\
                                    which have been remove adapter and primers. get all umi seqs \
                                    which length equal 14 nt.')
    parse.add_argument('-f1','--fastq1',required=True,help='forward fastq(R1)')
    parse.add_argument('-f2','--fastq2',required=True,help='reverse fastq(R2)')
    parse.add_argument('-o','--output',default='L_UMI_INFO.txt',help='output txt file. default=L_UMI_INFO.txt')
    args = parse.parse_args()
    
    logging.info(' L UMI fastq:R1=({}) R2=({})'.format(args.fastq1,args.fastq2))
    LR1_tmp=SeqIO.to_dict(SeqIO.parse(args.fastq1,'fastq'))
    LR2_tmp=SeqIO.to_dict(SeqIO.parse(args.fastq2,'fastq'))
    df_out = get_L_UMI(LR1_tmp,LR2_tmp)
    df_out.to_csv(args.output,sep='\t')