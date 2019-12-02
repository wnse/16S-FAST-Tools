#!/usr/bin/env python
# coding: utf-8
#from IPython.core.interactiveshell import InteractiveShell
#InteractiveShell.ast_node_interactivity = "all" 


import argparse
import logging
import os
import time
import sys
import re
import subprocess
import pysam
import shutil
import gzip

import random
import itertools
import multiprocessing

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


__version__ = "0.3.1_20191009"
parser = argparse.ArgumentParser()
parser.add_argument('--version', action='version',
                    version='%(prog)s {version}'.format(version=__version__))
parser.add_argument('-ar1','--Assemble_Read_1',help='Assemble Read 1 fastq data',required=True)
parser.add_argument('-ar2','--Assemble_Read_2',help='Assemble Read 2 fastq data',required=True)
parser.add_argument('-lr1','--Linker_Read_1',help='Linker Read 1 fastq data',required=True)
parser.add_argument('-lr2','--Linker_Read_2',help='Linker Read 2 fastq data',required=True)
parser.add_argument('-o','--out_dir',help='Result dir',required=True)
parser.add_argument('-n','--name',help='File name prefix',default='test')
parser.add_argument('-c','--CPUs',help='Number of CPUs. default 4',type=int,default=4)
parser.add_argument('-pc','--umi_pair_cutoff',\
                    help='Minimum number of paired umis appeared in Linker data. default 1',\
                    type=int,default=1)
parser.add_argument('-up','--unpaired_umi_seq',\
                    help='Group unpaired umi seq. 0:N 1:Y default 1',choices=[0,1],\
                    type=int,default=0)
parser.add_argument('-urc','--unpaired_rc_cutoff',\
                    help='Minimum unpaired umi reads count, used when unpaired_umi_seq is 0. default 200',\
                    type=int,default=200)
args = parser.parse_args()

umi_paire_min_counts = args.umi_pair_cutoff
unpaired_rc_cutoff = args.unpaired_rc_cutoff
umi_paire_less_cutoff = 0.5
con_seq_len_cutoff = 1200
threads = args.CPUs


# gunzip A1.fastq.gz
'''
A1_file = args.Assemble_Read_1
A1_fastq = ''
if os.path.splitext(A1_file)[1] == '.fastq':
    A1_fastq = A1_file
elif os.path.splitext(A1_file)[1] == '.fq':
    A1_fastq = A1_file
elif os.path.splitext(A1_file)[1] == '.gz':
    tmp_sys = os.system('gunzip ' + A1_file)
    A1_fastq = os.path.splitext(A1_file)[0]
if A1_fastq == '':
    sys.exit('A1_fastq is null')
'''

A1 = args.Assemble_Read_1
A2 = args.Assemble_Read_2
L1 = args.Linker_Read_1
L2 = args.Linker_Read_2

r1_adapter = 'TAGATCGC'
r2_adapter = 'CTAGTACG'
r1_primer = 'ATGGATGAGTCTGGGTGGAG'
r2_primer = 'ATCTTCATCTTTGCCCCCCT'

a_adapter = 'file:/Bioinfo/Database/16S_BSI/adpter/adapter.fa'
a_linker = 'file:/Bioinfo/Database/16S_BSI/adpter/linker.fa'
file_primer_rc = 'file:/Bioinfo/Database/16S_BSI/adpter/primer_rc.fa'
file_primer = 'file:/Bioinfo/Database/16S_BSI/adpter/primer.fa'
bowtie_db = '/Bioinfo/Database/16S_BSI/GM_BSI_v3'
taxon_map_file = '/Bioinfo/Database/16S_BSI/GM_BSI_v3_taxon_map.txt'
spades = '/root/anaconda3/bin/spades.py '
cutadapt = '/root/anaconda3/bin/cutadapt'
bowtie2 = '/root/anaconda3/bin/bowtie2'

def cutadapt_cmd(in_file,out_file,log_file,adapter,g_a):
    '''
    cutadapt remove fastq primer
    '''
    str_join=' '
    cmd=str_join.join([cutadapt,'-%s'%g_a,adapter,'-o',out_file,in_file, '-j', str(threads), '>',log_file])
    status=os.system(cmd)
    if status==0:
        logging.info(' done %s'%cmd)
    else:
        logging.info(' exit %s'%cmd)

def cutadapt_cmd_single(in_file,out_file,log_file,adapter,g_a):
    '''
    cutadapt remove fastq primer
    '''
    str_join=' '
    info_file = out_file + '_info_file'
    cmd=str_join.join([cutadapt,'-%s'%g_a,adapter,'-o',out_file,in_file,'--info-file',info_file,'>',log_file])
    status=os.system(cmd)
    if status==0:
        logging.info(' done %s'%cmd)
    else:
        logging.info(' exit %s'%cmd)

        
def bowtie_cmd(in_file,out_file,log_file):
    '''
    bowtie alignment
    '''
    str_join=' '
    cmd=str_join.join([bowtie2,'-p', str(threads),'-x',bowtie_db,'-f',in_file,'-S',out_file,'1>',log_file,'2>',log_file])
    status=os.system(cmd)
    if status==0:
        logging.info(' done %s'%cmd)
    else:
        logging.info(' exit %s'%cmd)

def reverse_complement(seq_list):
    '''
    sequences list reverse_complement
    '''
    rc_list=[]
    for seq in seq_list:
        rc_list.append(str(Seq(seq,generic_dna).reverse_complement()))
    return rc_list


def run_cmd(cmd,logfile):
    '''
    not useful any more
    '''
    with open(logfile,'w') as LOG:
        cmd_sub=subprocess.Popen(cmd,shell=True,stderr=LOG,stdout=LOG)
    return cmd_sub

def refresh_process(process_list):
    '''
    multiply process
    not useful any more
    '''
    aval=[]
    for i in process_list:
        if eval(i) == 0:
            #print ('{:>6}'.format(str(eval(i))),end='')
            aval.append(i)
        else:
            #print ('{:>6}'.format(str(eval(i).poll())),end='')
            aval.append(eval(i).poll())
    #print(' ===== ',end='')
    return aval

def Count_Circle(UMI_dict,UMI1_dict,UMI2_dict,cutoff=0.5,min_counts=10):
    ''' 
    1.排序UMI_dict中counts最多的idx；
    2.在剩余的idx中删除：
        已经确定的的UMI对的UMI2对应的所有其他UMI对中，该UMI2为最高数量的UMI1,除非该配对counts大于最多counts的cut_off倍数；
        已经确定的的UMI对的UMI1对应的所有其他UMI对中，该UMI1为最高数量的UMI2,除非该配对counts大于最多counts的cut_off倍数；
    3.如果配对counts小于min_counts，结束；
    4.如果无配对，递归结束；
    '''  
    max_tmp={}
    out_list=[]
    UMI_list = sorted(UMI_dict.items(),key=lambda UMI_dict:UMI_dict[1],reverse=True)
    for idx_count in UMI_list:
        if len(UMI1_dict) == 0:
            break
        (u1,u2)=idx_count[0].split('_')
        if u1 in UMI1_dict.keys():
            if len(UMI1_dict[u1]) == 0:
                del UMI1_dict[u1]
                continue
            if u2 not in UMI2_dict.keys():
                del UMI1_dict[u1][u2]
                continue
            if UMI1_dict[u1][u2] < min_counts:
                logging.info(' exit for {} counts {} < {}' \
                             .format(idx_count[0],UMI1_dict[u1][u2],min_counts))
                break
            out_list.append([idx_count[0],u1,u2,UMI1_dict[u1][u2]])
            list1=list(UMI2_dict[u2].keys())
            list2=list(UMI1_dict[u1].keys())
            max_counts=UMI1_dict[u1][u2]
            del UMI1_dict[u1]
            del UMI2_dict[u2]
            for deu in list1:
                if deu in UMI1_dict.keys():
                    if max(UMI1_dict[deu],key = UMI1_dict[deu].get) == u2:
                        if UMI1_dict[deu][u2] >= cutoff * max_counts:
                            out_list.append([idx_count[0],deu,u2,UMI1_dict[deu][u2]])                    
                        del UMI1_dict[deu]
            for deu in list2:
                if deu in UMI2_dict.keys():
                    if max(UMI2_dict[deu],key = UMI2_dict[deu].get) == u1:
                        if UMI2_dict[deu][u1] >= cutoff * max_counts:
                            out_list.append([idx_count[0],u1,deu,UMI2_dict[deu][u1]]) 
                        del UMI2_dict[deu]
    return out_list

def get_reads_counts(umi_list,counts_dict):
    count_list=[]
    for i in umi_list:
        if i in counts_dict.keys():
            count_list.append(counts_dict[i])
        else:
            count_list.append(0)
    return count_list

def iterator_slice(iterator, length):
    iterator = iter(iterator)
    while True:
        res = tuple(itertools.islice(iterator, length))
        if not res:
            break
        yield res 

def group_into_umis(seq):
    if seq.id in aUMI.keys():
        return (aUMI[seq.id],seq)
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


def queu_group_umi_seq(seq,umi2ID_dict,out_dir,paired=1):
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
            result = eval(func)(seq_tmp)
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
                        mkdir(dir_tmp_path)
                        umi_seq_dir_tmp_sub[umiID] = umi_seq_dir_tmp_sub.get(umiID,dir_tmp)
                        file_name =os.path.join(out_dir,umi_seq_dir_tmp_sub[umiID],umiID + '.fastq')
                        umi_file_list_sub[umi_seq_dir_tmp_sub[umiID]] = umi_file_list_sub.get(umi_seq_dir_tmp_sub[umiID],set())
                        umi_file_list_sub[umi_seq_dir_tmp_sub[umiID]].add(file_name)
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
                mkdir(dir_tmp_path)
                umi_seq_dir_tmp_sub[umiID] = umi_seq_dir_tmp_sub.get(umiID,dir_tmp)
                file_name =os.path.join(out_dir,umi_seq_dir_tmp_sub[umiID],umiID + '.fastq')
                umi_file_list_sub[umi_seq_dir_tmp_sub[umiID]] = umi_file_list_sub.get(umi_seq_dir_tmp_sub[umiID],set())
                umi_file_list_sub[umi_seq_dir_tmp_sub[umiID]].add(file_name)
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
    return (umi_counts_sub,umi_id_counts_sub,umi_file_list_sub,umi2seq_tmp) 
    

def spades_sub(file_list):
    file = file_list[0]
    outdir = file_list[1]
    logfile = file_list[2]
    cmd = spades + ' --s1 ' + file + ' -o ' + outdir + ' -t 1 ' \
                        +' 1>' + logfile + ' 2>' + logfile
    tmp_sys = os.system(cmd)
    rm_cmd = '/usr/bin/rm -rf outdir/assembly_graph* \
                            outdir/before_rr.fasta \
                            outdir/contigs.paths \
                            outdir/corrected \
                            outdir/dataset.info \
                            outdir/input_dataset.yaml \
                            outdir/misc \
                            outdir/K* \
                            outdir/params.txt \
                            outdir/tmp  \
                            outdir/scaffolds.* \
                            outdir/spades.log \
                            outdir/warnings.log \
                            outdir/configs '
    try:
        tmp_sys_rm = os.system(re.sub('outdir',outdir,rm_cmd))
    except OSError as e:
        logging.info('rm error {} - {}'.format(e.filename,e.strerror))
    return tmp_sys

def sam2table(sam_file):
    sam = pysam.AlignmentFile(sam_file, "r")
    sam_list=[]
    for i in sam:
        nm = 0
        if i.is_unmapped:
            sam_list.append([i.qname,i.qlen])
        else:
            for tmp in i.tags:
                if tmp[0] == 'NM':
                    nm = tmp[1]
            gap = len(i.blocks) - 1
            qlen = (i.qend - i.qstart + 1)
            iden = (qlen - gap - nm) / qlen
            sam_list.append([i.qname,i.qlen,i.qstart,i.qend,i.reference_name, \
                             i.reference_start,i.reference_end,i.mapq,nm,gap,iden,i.cigarstring])    
    df_sam = pd.DataFrame(sam_list,columns=(['qname','qlen','qstart','qend', 'rname','rstart', \
                                             'rend','mapq','nm','gap','identity','cigar']))
    df_taxon=pd.read_csv(taxon_map_file,sep='\t',names=['seq_id','taxon'],index_col=['seq_id'])
    df_total=pd.merge(df_sam,df_taxon,how='left',left_on='rname',right_index=True)
    df_total['taxon'].fillna('Unclassifid',inplace=True)
    return df_total

def assemble(file_list,out_dir,log_dir):
    '''
    组装序列
    调用 multiprocess.pool多进程
    使用map_async方法
    '''
    logging.info(' assemble {} start using threads:\t{}'.format(out_dir,threads))
    #tmp_file_list=list(file_list)[:]
    ass_success=0
    ass_failed=0
    for tmp_dir,file_name_list in file_list.items():
        queu_file=[]
        mkdir(os.path.join(out_dir,tmp_dir))
        mkdir(os.path.join(log_dir,tmp_dir))
        for f in file_name_list:
            f_name = os.path.splitext(os.path.basename(f))[0]
            spades_tmp_dir = os.path.join(out_dir,tmp_dir,f_name)
            logfile = os.path.join(log_dir,tmp_dir,f_name)
            queu_file.append([f,spades_tmp_dir,logfile])

        pool = multiprocessing.Pool(processes=threads)
        spades_return = pool.map_async(spades_sub,queu_file)
        spades_return.wait()

        if spades_return.ready():
            for get in spades_return.get():
                if get == 0:
                    ass_success += 1
                else:
                    ass_failed += 1
        logging.info(' {} end'.format(tmp_dir))
        
        pool.close()
        pool.join()

    logging.info(' assemble processed umi file sucessed {}:\t'.format(ass_success))
    logging.info(' assemble processed umi file failed {}:\t'.format(ass_failed))
    logging.info(' assemble end')

def merge_contigs(pass_file_name,filter_file_name,assemble_out_dir):
    # 处理组装后contig序列
    logging.info(' merge contig fasta start')
    umi_contig_count={}
    assemble_info={}
    contig_filtered=0
    contig_pass_filtered=0
    assemble_sucess=0
    assemble_failed=0
    handle_filter = open (filter_file_name,'a')
    handle_pass = open (pass_file_name,'a')
    for d in sorted(os.listdir(assemble_out_dir)):
        for d2 in os.listdir(os.path.join(assemble_out_dir,d)):
            tmp_path = os.path.join(assemble_out_dir,d,d2,'contigs.fasta')
            if os.path.exists(tmp_path): 
                tmp_fa = SeqIO.parse(tmp_path,'fasta')
                len_tmp = len(next(tmp_fa).seq)
                assemble_info[d2] = len_tmp
                assemble_sucess += 1
                if len_tmp < con_seq_len_cutoff:
                    tmp_fa = next(SeqIO.parse(tmp_path,'fasta'))
                    tmp_fa.id = d2 + '|' + tmp_fa.id
                    count_tmp = SeqIO.write(tmp_fa,handle_filter,'fasta')
                    umi_contig_count[d2] = umi_contig_count.get(d2,0)+count_tmp
                    contig_filtered += 1
                else:
                    tmp_filter_fa = []
                    tmp_fa = SeqIO.parse(tmp_path,'fasta')
                    for seq in tmp_fa:
                        if len(seq.seq) >= con_seq_len_cutoff:
                            contig_pass_filtered += 1
                            seq.id = d2 + '|' + seq.id
                            tmp_filter_fa.append(seq)
                    count_tmp=SeqIO.write(tmp_filter_fa,handle_pass,'fasta')
                    umi_contig_count[d2]=umi_contig_count.get(d2,0)+count_tmp
            else:
                assemble_info[d2] = 0
                assemble_failed += 1
        logging.info(' done merge contig of {}'.format(d))
    handle_filter.close()
    handle_pass.close()

    logging.info(' assemble successed:\t{}'.format(assemble_sucess))
    logging.info(' contig length pass filtered(>={}):\t{}'.format(con_seq_len_cutoff,contig_pass_filtered))
    logging.info(' contig length filtered(<{}):\t{}'.format(con_seq_len_cutoff,contig_filtered))
    #logging.info(' contig length filter Failed:')
    #logging.info(' {}'.format(contig_filtered))
    logging.info(' assemble failed:\t{}'.format(assemble_failed))
    #logging.info(assemble_failed)
    logging.info(' merge contig fasta end')
    return assemble_info,contig_pass_filtered
    
def sta_abundance(df_sub,identity,tag,pos):
    df_tmp = df_sub.loc[:]
    match="^" + tag[0] + "_"
    df_tmp[tag] = df_tmp['taxon'].str.split(';').str[pos].replace(regex={match:''})
    df_tmp.loc[df_tmp['taxon']=='Unclassifid',tag] = 'Unclassifid'
    df_tmp.loc[df_tmp['identity'] < identity,tag] = 'Unclassifid'
    abundance = pd.DataFrame(df_tmp.groupby(by=tag)['umi_id'].count().sort_values(ascending=False))
    abundance.columns=['Reads_Counts']
    abundance['Ratio%'] = abundance['Reads_Counts'] / abundance['Reads_Counts'].sum() *100
    return abundance

def format_taxon_file(df_sub,out_dir):
    ''' 
    file 
        taxon_abundance_all.csv 
        taxon_abundance_all.xlsx
    will write to out_dir
    '''
    df_tmp = df_sub.loc[:]
    df_tmp.sort_values(by='identity',inplace=True,ascending=False) # 按照identity由高到低排序
    df_tmp['umi_id'] = df_tmp['qname'].str.split('|',1).str[0] # 取出umi_id
    df_tmp.drop_duplicates(subset='umi_id',keep='first',inplace=True) # 相同umi_id只保留最前面的一个
    df_tmp = df_tmp[df_tmp['qlen']>=1200] # 保留序列长度大于1200bp的条目
    dict_tag = {'species':6,'genus':5,'family':4,'order':3,'class':2,'phylum':1,'kingdom':0} # 定义物种的位置
    iden_tag = {'species':0.97,'genus':0.945,'family':0.865,'order':0.82,'class':0.785,'phylum':0.75,'kingdom':0.6} # 定义物种的identity
    df_tmp.to_csv(os.path.join(out_dir,'taxon_alignment.csv'),sep='\t') # 筛选后的条目写入文件，用于提取fasta序列
    writer = pd.ExcelWriter(os.path.join(out_dir,'taxon_abundance.xlsx')) # 按照物种分类进行丰度统计，写入excel
    for tag,pos in dict_tag.items():
        df_tag = sta_abundance(df_tmp,iden_tag[tag],tag,pos)
        df_tag.to_excel(writer,tag)
    writer.save()
    return df_tmp['qname'].tolist()

def get_sub_fasta(sub_id_list,in_fa,out_fa):
    seq = SeqIO.to_dict(SeqIO.parse(in_fa,'fasta'))
    seq_sub = []
    for record in sub_id_list:
        seq_sub.append(seq[record])
    count_tmp = SeqIO.write(seq_sub,out_fa,'fasta')

def mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)
        logging.info(' mkdir '+ dir_name ) 

def plot_seq2UMI(np_seq2UMI,out_file):
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

def plot_contig_coverage(taxon_file,out_file):
    df = pd.read_csv(taxon_file, sep='\t',index_col=0)
    df.dropna(subset=['rname'],inplace=True)
    df['r_id'],df['r_len']=df['rname'].str.split('_').str
    df['r_len']=df['r_len'].apply(pd.to_numeric)
    df['s']=(1600 / df['r_len']) * df['rstart'] + 0.5
    df['e']=(1600 / df['r_len']) * df['rend'] + 0.5
    df.sort_values(by=['s','identity'],inplace=True)
    df.reset_index(inplace=True)
    array = np.zeros([df.shape[0],1600])
    for i in df.index:
        array[i][int(df.loc[i]['s']):int(df.loc[i]['e'])] = df.loc[i]['identity']   
    a4_dims = (11.7, 8.27)
    fig, axs = plt.subplots(2,figsize=a4_dims)
    sns.heatmap(array,ax=axs[1],cmap='Greys',xticklabels=100,yticklabels=10**(len(str(array.shape[0]))-1))
    sns.distplot(df['qlen'],ax=axs[0],bins=100,kde=False,color='black')
    axs[0].set(xlabel='Contig Length',ylabel='No. of Contigs')
    plt.tight_layout()
    plt.savefig(out_file,dpi=300)


def plot_contig_len(df_tmp,out_file):
    a4_dims = (11.7, 8.27)
    fig, axs = plt.subplots(figsize=a4_dims)
    df_tmp['log_reads_counts']=np.log10(df_tmp['reads_counts'])
    x=df_tmp[df_tmp['contig_len']>0]['log_reads_counts']
    y=df_tmp[df_tmp['contig_len']>0]['contig_len']
    axs.scatter(x=x,y=y,color='grey',alpha=0.8,s=0.5)
    # 拟合
    parameter = np.polyfit(x, y, 1)
    f = np.poly1d(parameter)
    axs.plot(x, f(x),color='black')
    x=df_tmp['log_reads_counts']
    y=df_tmp['contig_len']
    sns.lineplot(ax=axs,x=x,y=y,color='black',alpha=0.3)
    axs.axhline(1200,color='red',linestyle='-.',alpha=0.5)
    axs.set(xlabel='log10(Read Counts)', ylabel='Contig Length')
    plt.tight_layout()
    plt.savefig(out_file,dpi=600)

if not os.path.isdir(args.out_dir):
    os.makedirs(args.out_dir)
total_log = os.path.join(args.out_dir,'log')
logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO) #,filename = total_log,filemode = 'a')
logging.info(' version: {}'.format(__version__))
logging.info(' {}'.format(args.__dict__))

df_sta = pd.DataFrame(columns=['description','value'])
File_Tag=args.name
result_dir=os.path.join(args.out_dir,'result')
ana_dir=os.path.join(args.out_dir,'analysis')
mkdir(result_dir)
mkdir(ana_dir)

umi_seq_dir = os.path.join(ana_dir,'umi_seq')
spades_path = os.path.join(ana_dir,'spades_out')

A2_umi_file = os.path.join(ana_dir,File_Tag+'_cut_linker_P.fq')
spades_log_path = os.path.join(ana_dir,'spades_log')
umi1_fq = os.path.join(ana_dir,File_Tag+'_cut_primer_1.fq')
umi2_fq = os.path.join(ana_dir,File_Tag+'_cut_primer_2.fq')

con_seq_out_pass_raw_fa = os.path.join(ana_dir,'total_contig_pass.fasta')
con_seq_out_filter_raw_fa = os.path.join(ana_dir,'total_contig_filter.fasta')
con_seq_out_pass_fa = os.path.join(result_dir,'total_contig_pass.fasta')
con_seq_out_filter_fa = os.path.join(result_dir,'total_contig_filter.fasta')
pass_sam_file = os.path.join(ana_dir, 'total_contig_pass.bowtie.sam')
filter_sam_file = os.path.join(ana_dir, 'total_contig_filter.bowtie.sam')

# 去除连接文库adapter和linker，获得中间的UMI
#L1=os.path.join(rawdata_dir,'test_L_1.fq')
#L2=os.path.join(rawdata_dir,'test_L_2.fq')

logging.info(' cut L UMI start')
cutadapt_cmd(L1,os.path.join(ana_dir,File_Tag+'_cut_adapter_1.fq.gz'), \
             os.path.join(ana_dir,File_Tag+'_cut_adapter_1.log'),r1_adapter,'g')
cutadapt_cmd(L2,os.path.join(ana_dir,File_Tag+'_cut_adapter_2.fq.gz'), \
             os.path.join(ana_dir,File_Tag+'_cut_adapter_2.log'),r2_adapter,'g')

cutadapt_cmd(os.path.join(ana_dir,File_Tag+'_cut_adapter_1.fq.gz'), \
             umi1_fq,os.path.join(ana_dir,File_Tag+'_cut_primer_1.log'),r1_primer,'a')
cutadapt_cmd(os.path.join(ana_dir,File_Tag+'_cut_adapter_2.fq.gz'), \
             umi2_fq,os.path.join(ana_dir,File_Tag+'_cut_primer_2.log'),r2_primer,'a')


# 寻找拼接文库UMI

logging.info(' cut A UMI start')
cutadapt_cmd_single(A2,os.path.join(ana_dir,File_Tag+'_cut_adapter_P.fq.gz'), \
             os.path.join(ana_dir,File_Tag+'_cut_adapter_P.log'),a_adapter,'a')

cutadapt_cmd_single(os.path.join(ana_dir,File_Tag+'_cut_adapter_P.fq.gz'), \
             A2_umi_file,os.path.join(ana_dir,File_Tag+'_cut_linker_P.log'),a_linker,'g')

# 获得连接文库的R1和R2的UMI
logging.info(' UMI pairing start')
logging.info(' ' + umi1_fq + ' ' + umi2_fq)
LR1=SeqIO.to_dict(SeqIO.parse(umi1_fq,'fastq'))
LR2=SeqIO.to_dict(SeqIO.parse(umi2_fq,'fastq'))

umi1_Lreads={}
umi2_Lreads={}
umi1_LPreads={}
umi2_LPreads={}
pUMI1_dict={}
pUMI2_dict={}
pUMI_dict={}
umi_seq = 0
for idx in LR1.keys():
    seq1=str(LR1[idx].seq)
    seq2=str(LR2[idx].seq)
    u_idx = seq1 + '_' + seq2
    if (len(seq1) == 14) and (len(seq2) == 14):
        umi1_Lreads[seq1]=umi1_Lreads.get(seq1,0)+1
        umi2_Lreads[seq2]=umi2_Lreads.get(seq2,0)+1
        umi1_LPreads[seq1]=umi1_LPreads.get(seq1,0)+1
        umi2_LPreads[seq2]=umi2_LPreads.get(seq2,0)+1

        umi_seq += 1
        pUMI_dict[u_idx] = pUMI_dict.get(u_idx,0) + 1
        pUMI1_dict[seq1]=pUMI1_dict.get(seq1,{})
        pUMI2_dict[seq2]=pUMI2_dict.get(seq2,{})
        pUMI1_dict[seq1][seq2]=pUMI1_dict[seq1].get(seq2,0)+1  
        pUMI2_dict[seq2][seq1]=pUMI2_dict[seq2].get(seq1,0)+1 
    elif len(seq1) == 14:
        umi1_Lreads[seq1]=umi1_Lreads.get(seq1,0)+1
    elif len(seq2) == 14:
        umi2_Lreads[seq2]=umi2_Lreads.get(seq2,0)+1


# 获取UMI的配对关系
tmp1 = len(pUMI1_dict)
tmp2 = len(pUMI2_dict)
tmp3 = len(LR1)

logging.info(' No. of L seqs:\t{}'.format(tmp3))
logging.info(' No. of L seqs of UMIs:\t{}'.format(umi_seq))
logging.info(' No. of UMIs in L data(R1):\t{}'.format(tmp1))
logging.info(' No. of UMIs in L data(R2):\t{}'.format(tmp2))
df_sta = df_sta.append([{'description':'No. of L seqs','value':tmp3}],ignore_index=True)
df_sta = df_sta.append([{'description':'No. of L seqs of UMIs','value':umi_seq}],ignore_index=True)
df_sta = df_sta.append([{'description':'Ratio of L seqs of UMIs','value':umi_seq/tmp3}],ignore_index=True)
df_sta = df_sta.append([{'description':'No. of UMIs in L data','value':tmp1+tmp2}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tNo. of UMIs in L data(R1)','value':tmp1}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tRatio of UMIs in L data(R1)','value':tmp1/(tmp1+tmp2)}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tNo. of UMIs in L data(R2)','value':tmp2}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tRatio of UMIs in L data(R2)','value':tmp2/(tmp1+tmp2)}],ignore_index=True)

LR1=[]
LR2=[]
tmp_list=[]
tmp_list=Count_Circle(pUMI_dict,pUMI1_dict,pUMI2_dict,cutoff=umi_paire_less_cutoff, \
                      min_counts=umi_paire_min_counts)

pUMI_dict={}
pUMI1_dict={}
pUMI2_dict={}
df_umi_info = pd.DataFrame(tmp_list,columns=(['umi_id','umi1','umi2','pair_counts']))

for i in (['umi1_Lreads','umi1_LPreads']):
    df_tmp = pd.DataFrame.from_dict(eval(i),orient='index',columns=([i]))
    df_umi_info = pd.merge(df_umi_info,df_tmp,left_on='umi1',right_index=True,how='left')
umi1_Lreads={}
umi1_LPreads={}
for i in (['umi2_Lreads','umi2_LPreads']):
    df_tmp = pd.DataFrame.from_dict(eval(i),orient='index',columns=([i]))
    df_umi_info = pd.merge(df_umi_info,df_tmp,left_on='umi2',right_index=True,how='left')
umi2_Lreads={}
umi2_LPreads={}

df_umi_info.to_csv(os.path.join(result_dir,'umi_info.csv'),sep='\t')

tmp_list=[]

tmp1 = df_umi_info['umi_id'].unique().shape[0]
tmp2 = df_umi_info.shape[0]
tmp3 = df_umi_info['umi1'].unique().shape[0]
tmp4 = df_umi_info['umi2'].unique().shape[0]
logging.info(' No. of UMI ID(paired) in L data:\t{}({})'.format(tmp1,tmp2))
logging.info(' No. of paired UMIs in L data(R1):\t{}'.format(tmp3))
logging.info(' No. of paired UMIs in L data(R2):\t{}'.format(tmp4))
df_sta = df_sta.append([{'description':'No. of UMI ID(paired) in L data','value':tmp1,}],ignore_index=True)
df_sta = df_sta.append([{'description':'No. of paired UMIs in L data','value':tmp3+tmp4}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tNo. of paired UMIs in L data(R1)','value':tmp3}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tRatio of paired UMIs in L data(R1)','value':tmp3/(tmp3+tmp4)}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tNo. of paired UMIs in L data(R2)','value':tmp4}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tRatio of paired UMIs in L data(R2)','value':tmp4/(tmp3+tmp4)}],ignore_index=True)

df_umi_info['umi1_rc']=reverse_complement(df_umi_info['umi1'].tolist())
df_umi_info['umi2_rc']=reverse_complement(df_umi_info['umi2'].tolist())
umi2ID={**df_umi_info[['umi_id','umi1_rc']].set_index('umi1_rc').to_dict()['umi_id'], \
        **df_umi_info[['umi_id','umi2_rc']].set_index('umi2_rc').to_dict()['umi_id']}
logging.info(' UMI pairing end')

# 获取拼接文库中seq_id与umi对应关系
logging.info(' load A UMI data start')
aUMI={}
aUMI_unpaired={}
umi_seq=set()
umi_seq_unpaired={}
total_seq=0

total_seq = 0
n = 0
seq_type={}
for i in open(os.path.join(ana_dir,File_Tag+'_cut_adapter_P.fq.gz_info_file'),'r'):
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
df_sta = df_sta.append([{'description':'No. of A seq','value':total_seq}],ignore_index=True)

n = 0
m = 0 
for i in open(A2_umi_file+'_info_file','r'):
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
                    if seq in umi2ID.keys():
                        umi_seq.add(seq)
                        aUMI[seq_id]=seq
                    else:
                        umi_seq_unpaired[seq]=umi_seq_unpaired.get(seq,0)+1
                        aUMI_unpaired[seq_id]=seq
                seq_type[seq_id] = seq_type[seq_id] + '|' + tmp[7]
            else:
                del seq_type[seq_id]

logging.info(' No. of A reads with linkers:\t{}'.format(n))
type_seq={}
for k,v in seq_type.items():
    type_seq[v] = type_seq.get(v,0)+1

for k,v in type_seq.items():
    logging.info(' No. of A reads with adapter|linker({}):\t{}'.format(k,v))
    df_sta = df_sta.append([{'description':'No. of A reads with adapter|linker({})'.format(k),'value':v}],ignore_index=True)

# for AR in SeqIO.parse(A2_umi_file,'fastq'):
#     seq=str(AR.seq)
#     total_seq+=1
#     if len(seq)==14:
#         if seq in umi2ID.keys():
#             umi_seq.add(seq)
#             aUMI[AR.id]=seq
#         else:
#             umi_seq_unpaired[seq]=umi_seq_unpaired.get(seq,0)+1
#             aUMI_unpaired[AR.id]=seq         
            
logging.info(' load A_UMI Data end')
# 找到unpaired umi中大于 unpaired_rc_cutoff 的umi

df_unpaired_umi = pd.DataFrame.from_dict(umi_seq_unpaired,orient='index').reset_index()
df_unpaired_umi.columns = (['umi1_rc','reads_counts'])
df_unpaired_umi['umi_id'] = reverse_complement(df_unpaired_umi['umi1_rc'].tolist())

tmp1 = len(umi_seq)
tmp2 = df_unpaired_umi.shape[0]
tmp3 = len(aUMI.keys())
tmp4 = df_unpaired_umi['reads_counts'].sum()
logging.info(' No. of paired UMIs in A data:\t{:d}'.format(tmp1))
logging.info(' No. of unpaired UMI in A data:\t{:d}'.format(tmp2))
#logging.info(' No. of A seq:\t{:d}'.format(total_seq))
logging.info(' No. of A seq of paired UMI:\t{:d}({:.2f}%)'.format(tmp3,tmp3*100/total_seq))
logging.info(' No. of A seq of unpaired UMI:\t{:d}({:.2f}%)'.format(tmp4,tmp4*100/total_seq))
df_sta = df_sta.append([{'description':'No. of paired UMIs in A data','value':tmp1}],ignore_index=True)
df_sta = df_sta.append([{'description':'No. of unpaired UMIs in A data','value':tmp2}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tNo. of A seq of paired UMI','value':tmp3}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tRatio of A seq of paired UMI','value':tmp3/total_seq}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tNo. of A seq of unpaired UMI','value':tmp4}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tRatio of A seq of unpaired UMI','value':tmp4/total_seq}],ignore_index=True)

df_unpaired_umi = df_unpaired_umi[df_unpaired_umi['reads_counts'] > unpaired_rc_cutoff]. \
                    sort_values(by='reads_counts',ascending=False).reset_index()
df_unpaired_umi.drop(columns=['index'])
df_unpaired_umi.to_csv(os.path.join(result_dir,'umi_unpaired_info.csv'),sep='\t')
umi2ID_unpaired = df_unpaired_umi[['umi_id','umi1_rc']].set_index('umi1_rc').to_dict()['umi_id']
#df_sta.to_csv(os.path.join(result_dir,'analysis_sta.csv'),sep='\t',index=False)

# 通过配对UMI将连接文库序列分组
# 获取umi与umi_id的对应关系   
mkdir(umi_seq_dir)
mkdir(spades_path)
mkdir(spades_log_path)

paired_dir = os.path.join(umi_seq_dir,'paired')
mkdir(paired_dir)

logging.info(' write to umi paired seq file start')
(umi_counts,umi_id_counts,umi_file_list,umi2seq) = \
            queu_group_umi_seq(A1,umi2ID,paired_dir)
logging.info(' write to umi paired seq file end')


df_umi_info=pd.merge(df_umi_info,pd.DataFrame.from_dict(umi_id_counts,orient='index',\
                            columns=['reads_counts']),left_on='umi_id',right_index=True,how='left')
df_umi_info['umi1_reads_counts']=get_reads_counts(df_umi_info['umi1_rc'].tolist(),umi_counts)
df_umi_info['umi2_reads_counts']=get_reads_counts(df_umi_info['umi2_rc'].tolist(),umi_counts)
df_umi_info['reads_counts'].fillna(0,inplace=True)
df_umi_info.to_csv(os.path.join(result_dir,'umi_info.csv'),sep='\t')
logging.info(' No. of UMI_ID appeared in A data:\t{}'.format(umi2seq[-1][-1]))
df_sta = df_sta.append([{'description':'No. of UMI_ID in A data','value':umi2seq[-1][-1]}],ignore_index=True)
plot_seq2UMI(np.array(umi2seq),os.path.join(result_dir,'umi2seq.pdf'))
df_umi_seq = pd.DataFrame(np.array(umi2seq),columns=(['umi_seq','seq','umi','umi_id']))
df_umi_seq.to_csv(os.path.join(result_dir,'umi2seq.csv'),sep='\t')


assemble_out_dir=os.path.join(spades_path,'paired')
assemble_log_dir=os.path.join(spades_log_path,'paired')
mkdir(assemble_out_dir)
mkdir(assemble_log_dir)

assemble(umi_file_list,assemble_out_dir,assemble_log_dir)
assemble_info,contig_pass_filtered = merge_contigs(con_seq_out_pass_raw_fa,con_seq_out_filter_raw_fa,assemble_out_dir)

df_contig_len = pd.DataFrame.from_dict(assemble_info,orient='index',columns=['contig_len'])
df_contig_len = pd.merge(df_contig_len,df_umi_info[['umi_id','reads_counts']].drop_duplicates().set_index('umi_id'),left_index=True,right_index=True)
df_contig_len.to_csv(os.path.join(result_dir,'umi_contig_length.csv'),sep='\t')
plot_contig_len(df_contig_len.sort_values(by='contig_len'),os.path.join(result_dir,'umi_contig_length.png'))

tmp1=df_contig_len.shape[0]
tmp2=df_contig_len[df_contig_len['contig_len']>0].shape[0]
tmp3=df_contig_len[df_contig_len['contig_len']>=con_seq_len_cutoff].shape[0]
df_sta = df_sta.append([{'description':'\tNo. of contigs assembled','value':tmp2}],ignore_index=True)
df_sta = df_sta.append([{'description':'\tRatio of contigs assembled','value':tmp2/tmp1}],ignore_index=True)
df_sta = df_sta.append([{'description':'\t\tNo. of contigs filtered(>=1200,withPrimer)','value':tmp3}],ignore_index=True)
df_sta = df_sta.append([{'description':'\t\tRatio of contigs filtered(>=1200,withPrimer)','value':tmp3/tmp1}],ignore_index=True)

# 去除contig primer
# 比对contig seq
logging.info(' cut contig UMI primer start and mapping')
cutadapt_cmd(con_seq_out_pass_raw_fa, \
             os.path.join(ana_dir, 'total_contig_pass_cut_primerRC.fasta'), \
             os.path.join(ana_dir, 'total_contig_pass_cut_primerRC.log'),file_primer_rc, 'a')
cutadapt_cmd(os.path.join(ana_dir, 'total_contig_pass_cut_primerRC.fasta'), \
             con_seq_out_pass_fa, \
             os.path.join(ana_dir, 'total_contig_pass_cut_primerRC_cut_primer.log'), \
             file_primer, 'g')
bowtie_cmd(con_seq_out_pass_fa, \
           pass_sam_file,os.path.join(ana_dir, 'total_contig_pass_cut_primerRC_cut_primer.bowtie.log'))


cutadapt_cmd(con_seq_out_filter_raw_fa, \
             os.path.join(ana_dir, 'total_contig_filter_cut_primerRC.fasta'), \
             os.path.join(ana_dir, 'total_contig_filter_cut_primerRC.log'),file_primer_rc, 'a')
cutadapt_cmd(os.path.join(ana_dir, 'total_contig_filter_cut_primerRC.fasta'), \
             con_seq_out_filter_fa, \
             os.path.join(ana_dir, 'total_contig_filter_cut_primerRC_cut_primer.log'), \
             file_primer, 'g')
bowtie_cmd(con_seq_out_filter_fa, \
           filter_sam_file,os.path.join(ana_dir, 'total_contig_filter_cut_primerRC_cut_primer.bowtie.log'))


df_pass_total = sam2table(pass_sam_file)
df_filter_total = sam2table(filter_sam_file)
df_pass_total.to_csv(os.path.join(result_dir,'total_pass_bowtie_taxon.csv'),sep='\t')
df_filter_total.to_csv(os.path.join(result_dir,'total_filter_bowtie_taxon.csv'),sep='\t')
plot_contig_coverage(os.path.join(result_dir,'total_filter_bowtie_taxon.csv'),os.path.join(result_dir,'umi_contig_filter_cover.png'))

sub_fa_id = format_taxon_file(df_pass_total,result_dir)
get_sub_fasta(sub_fa_id,con_seq_out_pass_fa,os.path.join(result_dir,'taxon_seq.fasta'))
tmp1=len(sub_fa_id)
tmp2=df_contig_len.shape[0]
df_sta = df_sta.append([{'description':'\t\t\tNo. of contigs filtered(>=1200,withoutPrimer)','value':tmp1}],ignore_index=True)
df_sta = df_sta.append([{'description':'\t\t\tRatio of contigs filtered(>=1200,withoutPrimer)','value':tmp1/tmp2}],ignore_index=True)
df_sta.to_csv(os.path.join(result_dir,'analysis_sta.csv'),sep='\t',index=False)

if args.unpaired_umi_seq == 1:  
    unpaired_dir = os.path.join(umi_seq_dir,'unpaired')
    mkdir(unpaired_dir)
    
    logging.info(' write to umi unpaired seq file start')
    (umi_counts,umi_id_counts,umi_file_list,umi2seq) = \
                queu_group_umi_seq(A1,umi2ID_unpaired,unpaired_dir,paired=0)
    logging.info(' write to umi unpaired seq file end')
    
    df_unpaired_umi=pd.merge(df_unpaired_umi,pd.DataFrame.from_dict(umi_id_counts,orient='index',\
                                columns=['id_reads_counts']),left_on='umi_id',right_index=True,how='left')
    df_unpaired_umi['umi1_reads_counts']=get_reads_counts(df_unpaired_umi['umi1_rc'].tolist(),umi_counts)
    df_unpaired_umi['id_reads_counts'].fillna(0,inplace=True)
    df_unpaired_umi.to_csv(os.path.join(result_dir,'umi_unpaired_info.csv'),sep='\t')
    
    assemble_out_dir=os.path.join(spades_path,'unpaired')
    assemble_log_dir=os.path.join(spades_log_path,'unpaired')
    mkdir(assemble_out_dir)
    mkdir(assemble_log_dir)
    assemble(umi_file_list,assemble_out_dir,assemble_log_dir)
    
    con_seq_out_unpaired_pass_raw_fa = os.path.join(ana_dir,'total_contig_unpaired_pass.fasta')
    con_seq_out_unpaired_filter_raw_fa = os.path.join(ana_dir,'total_contig_unpaired_filter.fasta')
    con_seq_out_unpaired_pass_fa = os.path.join(result_dir,'total_contig_unpaired_pass.fasta')
    con_seq_out_unpaired_filter_fa = os.path.join(result_dir,'total_contig_unpaired_filter.fasta')
    assemble_sucess,contig_pass_filtered = merge_contigs(con_seq_out_unpaired_pass_fa,con_seq_out_unpaired_filter_fa,assemble_out_dir)

    if os.path.exists(con_seq_out_unpaired_pass_raw_fa):
        cutadapt_cmd(con_seq_out_unpaired_pass_raw_fa, \
                     os.path.join(ana_dir, 'total_contig_unpaired_pass_cut_primerRC.fasta'), \
                     os.path.join(ana_dir, 'total_contig_unpaired_pass_cut_primerRC.log'),file_primer_rc, 'a')
        cutadapt_cmd(os.path.join(ana_dir, 'total_contig_unpaired_pass_cut_primerRC.fasta'), \
                     con_seq_out_unpaired_pass_fa, \
                     os.path.join(ana_dir, 'total_contig_unpaired_pass_cut_primerRC_cut_primer.log'), \
                     file_primer, 'g')
        pass_sam_file = os.path.join(ana_dir, 'total_contig_unpaired_pass.bowtie.sam')
        bowtie_cmd(con_seq_out_unpaired_pass_fa, \
                   pass_sam_file, \
                   os.path.join(ana_dir, 'total_contig_unpaired_pass_cut_primerRC_cut_primer.bowtie.log'))
        df_pass_total = sam2table(pass_sam_file)
        df_pass_total.to_csv(os.path.join(result_dir,'total_unpaired_pass_bowtie_taxon.csv'),sep='\t')
    if os.path.exists(con_seq_out_unpaired_filter_raw_fa):
        cutadapt_cmd(con_seq_out_unpaired_filter_raw_fa, \
                     os.path.join(ana_dir, 'total_contig_unpaired_filter_cut_primerRC.fasta'), \
                     os.path.join(ana_dir, 'total_contig_unpaired_filter_cut_primerRC.log'),file_primer_rc, 'a')
        cutadapt_cmd(os.path.join(ana_dir, 'total_contig_unpaired_filter_cut_primerRC.fasta'), \
                     con_seq_out_unpaired_filter_fa, \
                     os.path.join(ana_dir, 'total_contig_unpaired_filter_cut_primerRC_cut_primer.log'), \
                     file_primer, 'g')
        filter_sam_file = os.path.join(ana_dir, 'total_contig_unpaired_filter.bowtie.sam')
        bowtie_cmd(con_seq_out_unpaired_filter_fa, \
                   filter_sam_file, \
                   os.path.join(ana_dir, 'total_contig_unpaired_filter_cut_primerRC_cut_primer.bowtie.log'))
        df_filter_total = sam2table(filter_sam_file)
        df_filter_total.to_csv(os.path.join(result_dir,'total_unpaired_filter_bowtie_taxon.csv'),sep='\t')


logging.info(' rmdir {}'.format(ana_dir))
shutil.rmtree(ana_dir)
logging.info(' End')
