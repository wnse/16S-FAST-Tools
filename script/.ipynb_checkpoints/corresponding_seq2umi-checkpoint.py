#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 10:27:33 2019

@author: yk
"""

import pandas as pd
import argparse
import logging
import sys
from . import reverse_complement

# 获取拼接文库中seq_id与umi对应关系

def corresponding_seq2umi(adapter_info_file,primer_info_file,umi_dict):
    sta_list = []
    aUMI = {}
    umi_seq = set()
    umi_seq_unpaired = {}
    total_seq = 0
    seq_type = {}
    
    n = 0
    for i in open(adapter_info_file):
        total_seq += 1
        tmp=i.split('\t',)
        seq_id = tmp[0].split(' ')[0]
        if tmp[1] == str(-1):
            continue
        else:
            n += 1
            seq_type[seq_id] = tmp[7]
    
    logging.info(' No. of A seq:\t{:d}'.format(total_seq))
    logging.info(' No. of A seq with adapters:\t{}'.format(n))
    sta_list.append(['No. of A seq',total_seq])
    
    n = 0
    m = 0 
    for i in open(primer_info_file):
        tmp=i.split('\t')
        seq_id = tmp[0].split(' ')[0]
        if seq_id in seq_type.keys():
            if tmp[1] == str(-1):
                del seq_type[seq_id]
            else:
                n += 1
                seq = tmp[6]
                if len(seq) == 14:
                    if seq_type[seq_id] == tmp[7]:
                        m += 1
                        if seq in umi_dict.keys():
                            aUMI[seq_id]=seq
                            umi_seq.add(seq)
                        else:
                            umi_seq_unpaired[seq]=umi_seq_unpaired.get(seq,0)+1
                    seq_type[seq_id] = seq_type[seq_id] + '|' + tmp[7]
                else:
                    del seq_type[seq_id]
    
    logging.info(' No. of A reads with linkers:\t{}'.format(n))
    type_seq={}
    for k,v in seq_type.items():
        type_seq[v] = type_seq.get(v,0)+1
    
    for k,v in type_seq.items():
        logging.info(' No. of A reads with adapter|linker({}):\t{}'.format(k,v))
        sta_list.append(['No. of A reads with adapter|linker({' + k + '})',v])
    
    df_unpaired_umi = pd.DataFrame.from_dict(umi_seq_unpaired,orient='index').reset_index()
    try:
        df_unpaired_umi.columns = (['umi1_rc','reads_counts'])
        tmp2 = df_unpaired_umi.shape[0]
        tmp4 = df_unpaired_umi['reads_counts'].sum()
    except ValueError as e:
        print (e)
        tmp4 = 0
        tmp2 = 0
        
    tmp1 = len(umi_seq)
    tmp3 = len(aUMI.keys())    
    logging.info(' No. of paired UMIs in A data:\t{:d}'.format(tmp1))
    logging.info(' No. of unpaired UMI in A data:\t{:d}'.format(tmp2))
    logging.info(' No. of A seq of paired UMI:\t{:d}({:.2f}%)'.format(tmp3,tmp3*100/total_seq))
    logging.info(' No. of A seq of unpaired UMI:\t{:d}({:.2f}%)'.format(tmp4,tmp4*100/total_seq))
    

    sta_list.append(['No. of paired UMIs in A data',tmp1])
    sta_list.append(['No. of unpaired UMIs in A data',tmp2])
    sta_list.append(['\tNo. of A seq of paired UMI',tmp3])
    sta_list.append(['\tRatio of A seq of paired UMI',tmp3/total_seq])
    sta_list.append(['\tNo. of A seq of unpaired UMI',tmp4])
    sta_list.append(['\tRatio of A seq of unpaired UMI',tmp4/total_seq])
    
    return aUMI,sta_list
    

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
    parser = argparse.ArgumentParser(description='match umi_id and seqs in A library \
                                     based on L_UMI_ID.txt file and A_R2 cut_adapter(primer)_info_file')
    parser.add_argument('-ai','--adapter_info',required=True,help='A R2 cut adapter info file')
    parser.add_argument('-pi','--primer_info',required=True,help='A R2 cut primer info file')
    parser.add_argument('-u','--umi_id',required=True,help='L UMI ID file')
    args = parser.parse_args()
    
    umi_id_file = args.umi_id
    df_umi_info = pd.read_csv(umi_id_file,sep='\t')
    df_umi_info['umi1_rc']=reverse_complement.reverse_complement(df_umi_info['u1'].tolist())
    df_umi_info['umi2_rc']=reverse_complement.reverse_complement(df_umi_info['u2'].tolist())
    umi2ID={**df_umi_info[['umi_id','umi1_rc']].set_index('umi1_rc').to_dict()['umi_id'], \
            **df_umi_info[['umi_id','umi2_rc']].set_index('umi2_rc').to_dict()['umi_id']}
    pd.DataFrame.from_dict(umi2ID,orient='index').to_csv('umi2umiID.txt',sep='\t')
    seq2umi,sta_list = corresponding_seq2umi(args.adapter_info,args.primer_info,umi2ID)
    pd.DataFrame.from_dict(seq2umi,orient='index').to_csv('seq2umi.txt',sep='\t',index=False)
    print(pd.DataFrame(sta_list))
    
