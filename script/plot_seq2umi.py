#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:08:05 2019

@author: yk
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_seq2umi(np_seq2UMI,out_file):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.set_xlabel('No. of Sequences')
    ax1.set_ylabel('No. of UMI Sequences')
    ax2.set_ylabel('No. of UMIs')
    x=np_seq2UMI[:,1]
    ax1.plot(x,np_seq2UMI[:,0], 'g-',label='UMI Sequences')   # green, solid line
    ax2.plot(x,np_seq2UMI[:,2], 'b-',label='UMIs') # blue
    ax2.plot(x,np_seq2UMI[:,3], 'r-',label='Paired UMIs') # blue
    ax1.legend(loc=8)
    ax2.legend(loc=4)
    plt.tight_layout()
    plt.savefig(out_file,dpi=600)

if __name__ == '__main__':
    parse = argparse.ArgumentParser(description='plot umi2seq file')
    parse.add_argument('-i','--input',required=True,help='input file(umi2seq.csv)')
    parse.add_argument('-o','--out',default='umi2seq.plot.pdf',help='output plot file.default=umi2seq.plot.pdf')
    args = parse.parse_args()
    
    a = np.array(pd.read_csv(args.input,sep='\t',index_col=0))
    plot_seq2UMI(a,args.out)
    