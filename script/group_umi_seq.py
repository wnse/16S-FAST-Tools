#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 16:24:11 2019

@author: yk
"""
import argparse
import logging
import os
import sys
import pandas as pd
import multiprocessing
from Bio import SeqIO
sys.path.append('/Users/yk/work/Microbiology/16S-FAST/')
from script_split import mkdir

def group_into_umis(seq,aUMI_dict):
    if seq.id in aUMI_dict.keys():
        return (aUMI_dict[seq.id],seq)
    else:
        return 0
    
def group_into_umis_unpaired(seq):
    if seq.id in aUMI_unpaired.keys():
        if aUMI_unpaired[seq.id] in umi2ID_unpaired.keys():
            return (aUMI_unpaired[seq.id],seq)
        else:
            return 0
    else:
        return 0

def write_to_file(tmp_list):
    count = 0
    file_name = tmp_list[0]
    file_name_seq_dict = tmp_list[1]
    with open(file_name,'a') as handle:
        for u in file_name_seq_dict.keys():
            seq = file_name_seq_dict[u]
            for record in seq:
                record.description=record.id + ' ' + u
            count += SeqIO.write(seq,handle,'fastq') 
    return count

def queu_group_umi_seq(seq,umi2ID_dict,out_dir,aUMI_dict,threads=1,paired=1):
    func='group_into_umis'
    if paired == 0:
        func='group_into_umis_unpaired'
        
    umi_counts_sub={}
    umi_id_counts_sub={}
    umi_file_list_sub={}
    umi_seq_dir_tmp_sub={}
    total_check_seq=0

    file_open = 'open'
    if os.path.splitext(seq)[1] == '.gz':
        file_open = 'gzip.open'

    with eval(file_open)(seq,'rt') as seq_handle:
        seqs = SeqIO.parse(seq_handle,'fastq')
        total_umi_seq_dict={}
        queue = []
        check_umiID={}
        count = 0
        count_seq = 0
        umi2seq_tmp=[]
        for seq_tmp in seqs:
            result = eval(func)(seq_tmp,aUMI_dict)
            count_seq += 1
            if result:
                umiID = umi2ID_dict[result[0]]
                total_umi_seq_dict[result[0]] = total_umi_seq_dict.get(result[0],[])
                total_umi_seq_dict[result[0]].append(result[1])
                total_check_seq += 1
                umi_counts_sub[result[0]] = umi_counts_sub.get(result[0],0)+1
                umi_id_counts_sub[umiID] = umi_id_counts_sub.get(umiID,0)+1
                if (len(total_umi_seq_dict.keys()) > 5000 or total_check_seq >= 500000):
                    file_name_seq_dict={}
                    queue_write = []
                    for u,record in total_umi_seq_dict.items():
                        umiID  = umi2ID_dict[u]
                        check_umiID[umiID] = check_umiID.get(umiID,0)+1
                        dir_tmp = 'tmp_' + str(int(len(check_umiID.keys()) / 1000))
                        dir_tmp_path = os.path.join(out_dir,dir_tmp) 
                        mkdir.mkdir(dir_tmp_path)
                        umi_seq_dir_tmp_sub[umiID] = umi_seq_dir_tmp_sub.get(umiID,dir_tmp)
                        file_name =os.path.join(out_dir,umi_seq_dir_tmp_sub[umiID],umiID + '.fastq')
                        file_name_seq_dict[file_name] = file_name_seq_dict.get(file_name,{})
                        file_name_seq_dict[file_name][u] = record
                    for file_name in file_name_seq_dict.keys():
                        queue_write.append([file_name,file_name_seq_dict[file_name]])
                    total_umi_seq_dict={}
                    file_name_seq_dict={}
                    total_check_seq=0
                    pool_write = multiprocessing.Pool(processes = threads) 
                    result_write = pool_write.map(write_to_file,queue_write)
                    pool_write.close()
                    pool_write.join()
                    for get in result_write:
                        count += get
                    logging.info(' write {} / {} seqs of {} umis in {} umi_IDs to files'. \
                                 format(count,count_seq,len(umi_counts_sub.keys()),len(umi_id_counts_sub.keys())))
                    umi2seq_tmp.append([count,count_seq,len(umi_counts_sub.keys()),len(umi_id_counts_sub.keys())])
    if total_umi_seq_dict:
        file_name_seq_dict={}
        queue_write = []
        for u,record in total_umi_seq_dict.items():
                umiID  = umi2ID_dict[u]
                check_umiID[umiID] = check_umiID.get(umiID,0)+1
                dir_tmp = 'tmp_' + str(int(len(check_umiID.keys()) / 1000))
                dir_tmp_path = os.path.join(out_dir,dir_tmp) 
                mkdir.mkdir(dir_tmp_path)
                umi_seq_dir_tmp_sub[umiID] = umi_seq_dir_tmp_sub.get(umiID,dir_tmp)
                file_name =os.path.join(out_dir,umi_seq_dir_tmp_sub[umiID],umiID + '.fastq')
                file_name_seq_dict[file_name] = file_name_seq_dict.get(file_name,{})
                file_name_seq_dict[file_name][u] = record
        for file_name in file_name_seq_dict.keys():
            queue_write.append([file_name,file_name_seq_dict[file_name]])
        total_umi_seq_dict={}
        file_name_seq_dict={}
        pool_write = multiprocessing.Pool(processes = threads) 
        result_write = pool_write.map(write_to_file,queue_write)
        pool_write.close()
        pool_write.join()
        for get in result_write:
            count += get
        logging.info(' write {} / {} seqs of {} umis in {} umi_IDs to files'. \
                     format(count,count_seq,len(umi_counts_sub.keys()),len(umi_id_counts_sub.keys())))
        umi2seq_tmp.append([count,count_seq,len(umi_counts_sub.keys()),len(umi_id_counts_sub.keys())])
    return (umi_counts_sub,umi_id_counts_sub,umi2seq_tmp) 
    
if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
    parse = argparse.ArgumentParser(description='write each A seq to corresponded umi id file')
    parse.add_argument('-s','--seq2umi',help='seq id and umi table file(seq2umi.txt)',required=True)
    parse.add_argument('-u','--umi_id',required=True,help='L UMI ID file')
    parse.add_argument('-f','--fastq',help='A fastq file',required=True)
    parse.add_argument('-o','--outdir',help='output dir',required=True)
    parse.add_argument('-t','--threads',help='cpu threads,default=1',default=1,type=int)
    args = parse.parse_args()
    

    mkdir.mkdir(args.outdir)
    df_aumi = pd.read_csv(args.seq2umi,sep='\t')
    df_aumi.columns = (['seqid','umi'])
    aUMI = dict(zip(df_aumi['seqid'],df_aumi['umi']))
    umi_id_file = args.umi_id
    df_umi_info = pd.read_csv(umi_id_file,sep='\t')
    df_umi_info['u1_rc']=reverse_complement.reverse_complement(df_umi_info['u1'].tolist())
    df_umi_info['u2_rc']=reverse_complement.reverse_complement(df_umi_info['u2'].tolist())
    umi2ID={**df_umi_info[['umi_id','u1_rc']].set_index('u1_rc').to_dict()['umi_id'], \
            **df_umi_info[['umi_id','u2_rc']].set_index('u2_rc').to_dict()['umi_id']}
    
    (umi_counts,umi_id_counts,umi_file_list,umi2seq)= \
        queu_group_umi_seq(args.fastq,umi2ID,args.outdir,aUMI,args.threads)
    
    df_umi_info=pd.merge(df_umi_info,pd.DataFrame.from_dict(umi_id_counts,orient='index',\
                                columns=['umi_id_reads_counts']),left_on='umi_id',right_index=True,how='left')
    
    df_tmp = pd.DataFrame.from_dict(umi_counts,orient='index',columns=['reads_counts'])
    df_umi_info=pd.merge(df_umi_info,df_tmp,left_on='u1_rc',right_index=True,how='left')
    df_umi_info.rename(columns={'reads_counts':'umi1_reads_counts'},inplace=True)
    df_umi_info=pd.merge(df_umi_info,df_tmp,left_on='u2_rc',right_index=True,how='left')
    df_umi_info.rename(columns={'reads_counts':'umi2_reads_counts'},inplace=True)    
    
    df_umi_info.fillna(0,inplace=True)
    df_umi_info.to_csv(args.umi_id + '.ReadsCounts.txt',sep='\t')
        
    