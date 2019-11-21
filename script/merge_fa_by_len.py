#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 00:26:46 2019

@author: yk
"""
import argparse
import os
from Bio import SeqIO

def merge_fa_by_len(indir,out_pass,out_fail,cutoff,min_value=0):
    contigs = 0
    nocontigs =0
    pass_no = 0
    fail_no = 0
    min_no = 0
    pass_fa = []
    fail_fa = []
    s_dirs = os.listdir(indir)
    total_no = len(s_dirs)
    dir_path = map(lambda x:os.path.join(indir,x),s_dirs)
    for d in dir_path:
        file = os.path.join(d,'contigs.fasta')
        contigs += 1
        if os.path.exists(file):
            seqs = SeqIO.parse(file,'fasta')
            len_list = map(lambda x: len(x.seq),seqs)
            if max(len_list) >= cutoff:
                pass_no += 1
                for seq in SeqIO.parse(file,'fasta'):
                    n = 0
                    if len(seq.seq) >= cutoff:
                        seq.id = d + '_' + str(n)
                        pass_fa.append(seq)
                        n += 1
            elif max(len_list) <= min_value:
                min_no += 1
            else:
                fail_no += 1
                for seq in SeqIO.parse(file,'fasta'):
                    n = 0
                    if len(seq.seq) > min_value:
                        seq.id = d + '_' + str(n)
                        fail_fa.append(seq)
                        n += 1
        else:
            nocontigs += 1
    pass_seq = SeqIO.write(out_pass,pass_fa,'fasta')
    fail_seq = SeqIO.write(out_fail,fail_fa,'fasta')
    print('contigs:\t{}'.format(contigs))
    print('nocontigs:\t{}'.format(nocontigs))
    print('pass no:\t{}'.format(pass_no))
    print('pass seqs:\t{}'.format(pass_seq))
    print('fail no:\t{}'.format(fail_no))
    print('fail) seqs:\t{}'.format(fail_seq))
    print('min value seq:\t{}'.format(min_no))
            
if __name__ == '__main__':
    indir = '/Users/yk/work/Microbiology/16S-FAST/script_split/test/analysis/spades/tmp_0/'
    out_pass = 'pass.fa'
    out_fail = 'fail.fa'
    cutoff = 1200
    min_value = 100
    merge_fa_by_len(indir,out_pass,out_fail,cutoff,min_value)
               
            