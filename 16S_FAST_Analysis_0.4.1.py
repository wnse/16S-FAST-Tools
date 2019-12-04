#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:50:33 2019

@author: yk
"""

import argparse
import logging
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

from script import *

logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
__version__ = "0.3.1_20191009"
logging.info(' version: {}'.format(__version__))
logging.info(' script path:{}'.format(sys.path[0]))

parser = argparse.ArgumentParser()
parser.add_argument('--version', action='version',
                    version='%(prog)s {version}'.format(version=__version__))
parser.add_argument('-ar1','--Assemble_Read_1',required=True,
                    help='Assemble Read 1 fastq data')
parser.add_argument('-ar2','--Assemble_Read_2',required=True,
                    help='Assemble Read 2 fastq data')
parser.add_argument('-lr1','--Linker_Read_1',required=True,
                    help='Linker Read 1 fastq data')
parser.add_argument('-lr2','--Linker_Read_2',required=True,
                    help='Linker Read 2 fastq data')
parser.add_argument('-o','--out_dir',required=True,
                    help='Result dir')
parser.add_argument('-n','--name',default='test',
                    help='File name prefix',)
parser.add_argument('-c','--CPUs',type=int,default=4,
                    help='Number of CPUs. default 4',)
parser.add_argument('-pc','--umi_pair_cutoff',type=int,default=1,
                    help='Minimum number of paired umis appeared in Linker '
                    'data.',)
parser.add_argument('-up','--unpaired_umi_seq',type=int,default=0,choices=[0,1],
                    help='Group unpaired umi seq. 0:N 1:Y default 1',)
parser.add_argument('-urc','--unpaired_rc_cutoff',type=int,default=200,
                    help='Minimum unpaired umi reads count, used when '
                    'unpaired_umi_seq is 0.',)
parser.add_argument('-minl','--minlength',default=1200,type=int,
                    help='Minimum length of contigs filtered')
parser.add_argument('-maxl','--maxlength',default=0,type=int,
                    help='Maximum length of contigs filtered,0 means no limits')

raw_data_dir = '/root/test/rawdata/'
output_dir = '/root/test/test'
args = parser.parse_args(['-ar1',os.path.join(raw_data_dir,'FL-P_BKDL190814481-1a-AK3457-A26_1.fq.gz'),
                          '-ar2',os.path.join(raw_data_dir,'FL-P_BKDL190814481-1a-AK3457-A26_2.fq.gz'),
                          '-lr1',os.path.join(raw_data_dir,'FL-L_BKDL190814481-1a-AK2926-A26_1.fq.gz'),
                          '-lr2',os.path.join(raw_data_dir,'FL-L_BKDL190814481-1a-AK2926-A26_2.fq.gz'),
                          '-o',output_dir,
                          '-pc','50',
                          '-c','7',
                          '-minl','1200',
                          '-maxl','1700'])

logging.info(' {}'.format(args.__dict__))

umi_paire_min_counts = args.umi_pair_cutoff
unpaired_rc_cutoff = args.unpaired_rc_cutoff
umi_paire_less_cutoff = 0.5
con_seq_len_cutoff = 1200
threads = args.CPUs

A1 = args.Assemble_Read_1
A2 = args.Assemble_Read_2
L1 = args.Linker_Read_1
L2 = args.Linker_Read_2

cmd_path = '/Bioinfo/bin/'
bin_path = '/root/anaconda3/bin/'
db_path = '/Bioinfo/Database/'

bowtie_db = db_path + '16S_BSI/GM_BSI_v3'
bowtie_db_taxon_map_file = db_path + 'GM_BSI_v3_taxon_map.txt'
mothur_db_fa = db_path + 'mothur/silva_132_v3.fa'
mothur_db_tax = db_path + 'mothur/silva_132_v3.2.tax'

spades = bin_path + 'spades.py'
bowtie2 = bin_path + 'bowtie2'
cutadapt = bin_path + 'cutadapt'

mothur = cmd_path + 'mothur'

a_adapter = ('file:' + os.path.join(sys.path[0] + '/adapter/adapter.fa'))
a_linker = ('file:' + os.path.join(sys.path[0] + '/adapter/linker.fa'))
file_primer_rc = ('file:' + os.path.join(sys.path[0] + '/adapter/primer_rc.fa'))
file_primer = ('file:' + os.path.join(sys.path[0] + '/adapter/primer.fa'))
r1_adapter = 'TAGATCGC'
r2_adapter = 'CTAGTACG'
r1_primer = 'ATGGATGAGTCTGGGTGGAG'
r2_primer = 'ATCTTCATCTTTGCCCCCCT'
                       
sta_list = []
df_sta = pd.DataFrame(columns=['description','value'])
File_Tag = args.name
result_dir = os.path.join(args.out_dir,'result')
ana_dir = os.path.join(args.out_dir,'analysis')
mkdir.mkdir(result_dir)
mkdir.mkdir(ana_dir)

# 去除连接文库adapter和linker，获得中间的UMI
L1_cut_adapter = os.path.join(ana_dir,File_Tag+'_cut_adapter_1.fq.gz')
L1_cut_adapter_log = os.path.join(ana_dir,File_Tag+'_cut_adapter_1.log')
L2_cut_adapter = os.path.join(ana_dir,File_Tag+'_cut_adapter_2.fq.gz')
L2_cut_adapter_log = os.path.join(ana_dir,File_Tag+'_cut_adapter_2.log')

umi1_fq = os.path.join(ana_dir,File_Tag+'_cut_primer_1.fq')
umi1_fq_log = os.path.join(ana_dir,File_Tag+'_cut_primer_1.log')
umi2_fq = os.path.join(ana_dir,File_Tag+'_cut_primer_2.fq')
umi2_fq_log = os.path.join(ana_dir,File_Tag+'_cut_primer_2.log')

logging.info(' cut L UMI start')
submit_cutadapt.submit_cutadapt(L1,L1_cut_adapter,L1_cut_adapter_log,
                                r1_adapter,'g',cutadapt,threads)
submit_cutadapt.submit_cutadapt(L2,L2_cut_adapter,L2_cut_adapter_log,
                                r2_adapter,'g',cutadapt,threads)

submit_cutadapt.submit_cutadapt(L1_cut_adapter,umi1_fq,umi1_fq_log,
                                r1_primer,'a',cutadapt,threads)
submit_cutadapt.submit_cutadapt(L2_cut_adapter,umi2_fq,umi2_fq_log,
                                r2_primer,'a',cutadapt,threads)

# 寻找拼接文库UMI
A2_cut_adapter = os.path.join(ana_dir,File_Tag+'_cut_adapter_P.fq.gz')
A2_cut_adapter_log = os.path.join(ana_dir,File_Tag+'_cut_adapter_P.log')
A2_cut_adapter_info = A2_cut_adapter + '.cutadapt.info.file'
A2_umi_file = os.path.join(ana_dir,File_Tag+'_cut_linker_P.fq')
A2_umi_file_log = os.path.join(ana_dir,File_Tag+'_cut_linker_P.log')
A2_umi_file_info = A2_umi_file + '.cutadapt.info.file'

logging.info(' cut A UMI start')
submit_cutadapt.submit_cutadapt(A2,A2_cut_adapter,A2_cut_adapter_log,
                                a_adapter,'a',cutadapt,info=1)
submit_cutadapt.submit_cutadapt(A2_cut_adapter,A2_umi_file,A2_umi_file_log,
                                a_linker,'g',cutadapt,info=1)

# 获得连接文库的R1和R2的UMI
logging.info(' UMI pairing start')

LR1=SeqIO.to_dict(SeqIO.parse(umi1_fq,'fastq'))
LR2=SeqIO.to_dict(SeqIO.parse(umi2_fq,'fastq'))
L_UMI_INFO_FILE = os.path.join(result_dir,File_Tag+'.L.umi.INFO.txt')
L_UMI_ID_FILE = os.path.join(result_dir,File_Tag+'.L.umiID.INFO.txt')
SEQ2UMI_FILE = os.path.join(result_dir,File_Tag+'.seq2umi.txt')

df_L_UMI_INFO = get_L_UMI.get_L_UMI(LR1,LR2)
LR1 = {}
LR2 = {}
df_L_UMI_INFO.to_csv(L_UMI_INFO_FILE,sep='\t')
df_L_UMI_INFO = df_L_UMI_INFO[~ ((df_L_UMI_INFO['u1'].str.contains('\*')) 
                                 |(df_L_UMI_INFO['u2'].str.contains('\*')))]
umi = {}
umi1 = {}
umi2 = {}
for i in df_L_UMI_INFO.index.to_list():
    u1 = str(df_L_UMI_INFO.loc[i,'u1'])
    u2 = str(df_L_UMI_INFO.loc[i,'u2'])
    u = u1 + '_' + u2
    c = int(df_L_UMI_INFO.loc[i,'counts_of_paire'])
    umi[u] = c
    umi1[u1] = umi1.get(u1,{})
    umi1[u1][u2] = c
    umi2[u2] = umi2.get(u2,{})
    umi2[u2][u1] = c

tmp_list = define_UMI_ID.define_UMI_ID(umi,umi1,umi2,
                                       min_counts=args.umi_pair_cutoff)
if tmp_list:
    df_L_UMI_INFO = pd.DataFrame(tmp_list)
    df_L_UMI_INFO.columns=(['umi_id','u1','u2','counts_of_paire'])
    df_L_UMI_INFO.to_csv(L_UMI_ID_FILE,sep='\t',index=False)
    df_L_UMI_INFO['u1_rc']=reverse_complement.reverse_complement(
            df_L_UMI_INFO['u1'].tolist())
    df_L_UMI_INFO['u2_rc'] = reverse_complement.reverse_complement(
            df_L_UMI_INFO['u2'].tolist())
    umi2ID={
        **df_L_UMI_INFO[['umi_id','u1_rc']].set_index('u1_rc').to_dict()['umi_id'],
        **df_L_UMI_INFO[['umi_id','u2_rc']].set_index('u2_rc').to_dict()['umi_id']
        }
else:
    sys.exit(' No UMI ID defined !')
    
seq2umi,tmp_list = corresponding_seq2umi.corresponding_seq2umi(
    A2_cut_adapter_info,
    A2_umi_file_info,
    umi2ID)
sta_list.extend(tmp_list)
pd.DataFrame.from_dict(seq2umi,orient='index').to_csv(SEQ2UMI_FILE,sep='\t')

# 分组A R1序列
umi_seq_dir = os.path.join(ana_dir,'umi_seq')
mkdir.mkdir(umi_seq_dir)
logging.info(' write to umi paired seq file start')
(umi_counts,umi_id_counts,umi2seq) = group_umi_seq.queu_group_umi_seq(
        A1,umi2ID,umi_seq_dir,seq2umi,threads)
logging.info(' write to umi paired seq file end')

df_L_UMI_INFO = pd.merge(
        df_L_UMI_INFO,
        pd.DataFrame.from_dict(umi_id_counts,
                               orient = 'index',
                               columns = ['umi_id_reads_counts']),
        left_on='umi_id',
        right_index=True,
        how='left')

df_tmp = pd.DataFrame.from_dict(
        umi_counts,
        orient='index',
        columns=['reads_counts'])
df_L_UMI_INFO = pd.merge(
        df_L_UMI_INFO,
        df_tmp,
        left_on='u1_rc',
        right_index=True,
        how='left')
df_L_UMI_INFO.rename(
        columns={'reads_counts':'umi1_reads_counts'},
        inplace=True)
df_L_UMI_INFO = pd.merge(
        df_L_UMI_INFO,
        df_tmp,
        left_on='u2_rc',
        right_index=True,how='left')
df_L_UMI_INFO.rename(
        columns={'reads_counts':'umi2_reads_counts'},
        inplace=True)
df_L_UMI_INFO.fillna(0,inplace=True) 
df_L_UMI_INFO.to_csv(L_UMI_ID_FILE,sep='\t',index=False)

logging.info(' No. of UMI_ID appeared in A data:\t{}'.format(umi2seq[-1][-1]))
sta_list.append(['No.of UMI_ID appeared in A data:',umi2seq[-1][-1]])

umi2seq_fig = os.path.join(result_dir,File_Tag+'.umi2seq.pdf')
umi2seq_csv = os.path.join(result_dir,File_Tag+'.umi2seq.csv')
plot_seq2umi.plot_seq2umi(np.array(umi2seq),umi2seq_fig)
pd.DataFrame(np.array(umi2seq),
             columns=(['umi_seq','seq','umi','umi_id'])).to_csv(
                     umi2seq_csv,sep='\t')

umi_id_fastq_file_list = os.path.join(result_dir,File_Tag+'.umi.fastq.list.txt')

# 组装序列
spades_path = os.path.join(ana_dir,'spades')
spades_log_path = os.path.join(ana_dir,'spades_log')
mkdir.mkdir(spades_path)
mkdir.mkdir(spades_log_path)
tmp_list=[]

for tmp_dir in os.listdir(umi_seq_dir):
    spades_path_tmp = os.path.join(spades_path,tmp_dir)
    spades_log_tmp = os.path.join(spades_log_path,tmp_dir)
    mkdir.mkdir(spades_path_tmp)
    mkdir.mkdir(spades_log_tmp)
    file_path_tmp = os.path.join(umi_seq_dir,tmp_dir)
    files = map(lambda x:os.path.join(file_path_tmp,x),
                os.listdir(file_path_tmp))
    logging.info(' assemble {} files from {} using {} by {} threads:\t'
                 .format(len(os.listdir(file_path_tmp)),
                         file_path_tmp,spades,threads))
    tmp_list.append(
            assemble_fastq_from_list.assemble(
                    files,
                    spades_path_tmp,
                    spades_log_tmp,
                    spades,threads))

sta_list.append(['No.of UMI_ID success assembled:',tmp_list[0]])

merge_fa = os.path.join(result_dir,'merged.fasta')
merge_fa_info = os.path.join(result_dir,'merged.fasta.info')
merge_filter_fa = os.path.join(result_dir,'merged.filter.fasta')
ID_info = os.path.join(result_dir ,'asv.id.info')
asv_fa = os.path.join(result_dir,'asv.fasta')
assemble_list = []
for tmpdir in os.listdir(spades_path):
    for d in os.listdir(os.path.join(spades_path, tmpdir)):
        contig_file = os.path.join(spades_path, tmpdir, d, 'contigs.fasta')
        if os.path.exists(contig_file):
            tmp = [d,contig_file]
            assemble_list.append(tmp)
    
df_merge_fa = merge_contigs.merge_contigs(assemble_list,merge_fa)
df_merge_fa.to_csv(merge_fa_info,sep='\t')
tmp, df_cut_info = cut_fa_by_len.cut_fa_by_len(merge_fa,merge_filter_fa,
                                               args.minlength,args.maxlength)
sta_list.append(['No. Contigs filtered(1200-1700bp):',tmp])
df_cut_info.to_csv(ID_info,sep='\t')

submit_mothur.submit_mothur(merge_filter_fa,mothur_db_fa,mothur_db_tax,mothur,threads)
tax_file = (re.search(r'(.*)\.fasta',merge_filter_fa).group(1) + 
              re.search(r'(\..*)\.tax',mothur_db_tax).group(1) + 
              '.wang.taxonomy')
df_consensus, df_unconsensus = get_consensus_seq_from_mothur.get_consensus_seq_from_mothur(
    tax_file,merge_filter_fa)

'''
df_asv = get_asv_seq_from_fasta.get_asv_seq_from_fasta(
    merge_filter_fa,args.name,asv_fa,'asv')
df_asv = pd.merge(
    df_cut_info,
    df_asv,
    left_index=True,
    right_on='Seq_ID').set_index('Seq_ID')
df_asv.to_csv(ID_info,sep='\t')
'''
print(pd.DataFrame(sta_list))





