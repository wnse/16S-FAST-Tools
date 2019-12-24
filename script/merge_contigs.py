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


def merge_contigs(contig_list, merge_fa_name):
    '''合并组装后contigs
    参数：
        contig_list: contig文件路径列表（type=list）
        merge_fa_name: 输出文件
    返回：
        df: 每个umiID的contig数目（type=pandas）
    '''
    seq_no_dict = {}
    with open(merge_fa_name, 'w') as w:
        m = 0
        for record in contig_list:
            name = record[0]
            fa = record[1]
            tmp = []
            n = 0
            for seq in SeqIO.parse(fa, 'fasta'):
                seq.description = str(name) + "|" + str(seq.description)
                seq.id = str(m) + '_' + str(n)
                tmp.append(seq)
                n += 1
            seq_no = SeqIO.write(tmp, w, 'fasta')
            seq_no_dict[name] = [m, seq_no]
            m += 1
    return pd.DataFrame.from_dict(seq_no_dict, orient='index',
                                  columns=(['umi_id', 'NO_of_Contigs']))


if __name__ == '__main__':
    parse = argparse.ArgumentParser(description='merge contig fastas in a dir',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-d', '--indir', required=True,
                       help='input dir.input dir-dir-fasta')
    parse.add_argument('-o', '--outfasta', default='merged.fasta',
                       help='output merge.fasta')
    parse.add_argument('-i', '--info', default='merged.fasta.info',
                       help='output number of contigs info of each sample')
    args = parse.parse_args()

    input_list = []
    for d in os.listdir(args.indir):
        contig_file = os.path.join(args.indir, d, 'contigs.fasta')
        if os.path.exists(contig_file):
            tmp = [d, contig_file]
            input_list.append(tmp)

    df = merge_contigs(input_list, args.outfasta)
    df.to_csv(args.info, sep='\t')
