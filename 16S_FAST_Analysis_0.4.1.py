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
import shutil
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
from script import *

logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
__version__ = "0.4_20191205"
logging.info(' version: {}'.format(__version__))
logging.info(' script path:{}'.format(sys.path[0]))

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
                    help='Number of CPUs.',)
parser.add_argument('-pc','--umi_pair_cutoff',type=int,default=1,
                    help='Minimum number of paired umis appeared in Linker '
                    'data.',)
parser.add_argument('-minl','--minlength',default=1200,type=int,
                    help='Minimum length of contigs filtered')
parser.add_argument('-maxl','--maxlength',default=0,type=int,
                    help='Maximum length of contigs filtered,0 means no limits')
parser.add_argument('-group','--group',action='store_true',help='keep group umi seq in result')
parser.add_argument('-debug','--debug',action='store_true',help='not delete all processing files')
args = parser.parse_args()
#raw_data_dir = '/mnt/data/work/test_16S_FAST_v4/rawdata/'
#output_dir = '/mnt/data/work/test_16S_FAST_v4/test_out_1211'
#args = parser.parse_args(['-ar1',os.path.join(raw_data_dir,'test_P_1.fq.gz'),
#                          '-ar2',os.path.join(raw_data_dir,'test_P_2.fq.gz'),
#                          '-lr1',os.path.join(raw_data_dir,'test_L_1.fq.gz'),
#                          '-lr2',os.path.join(raw_data_dir,'test_L_2.fq.gz'),
# raw_data_dir = '/mnt/data/16S-FAST_data/20191209/rawdata/'  
# output_dir = '/mnt/data/work/test_16S_FAST_v4/test_out_LJ_new'
# args = parser.parse_args(['-ar1',os.path.join(raw_data_dir,'LJ-P_S33_L008_R1_001.fastq.gz'),
#                           '-ar2',os.path.join(raw_data_dir,'LJ-P_S33_L008_R2_001.fastq.gz'),
#                           '-lr1',os.path.join(raw_data_dir,'LJ-L_S34_L008_R1_001.fastq.gz'),
#                           '-lr2',os.path.join(raw_data_dir,'LJ-L_S34_L008_R2_001.fastq.gz'), 
#                           '-o',output_dir,
#                           '-pc','3',
#                           '-c','24',
#                           '-minl','800',
#                           '-maxl','1700',
#                           '-n','LJ'])

logging.info(' {}'.format(args.__dict__))

umi_paire_min_counts = args.umi_pair_cutoff
threads = args.CPUs

A1 = args.Assemble_Read_1
A2 = args.Assemble_Read_2
L1 = args.Linker_Read_1
L2 = args.Linker_Read_2

cmd_path = '/Bioinfo/bin/'
bin_path = '/root/anaconda3/bin/'
db_path = '/Bioinfo/Database/Silva_132_v3/'

bowtie_db = db_path + '16S_BSI/GM_BSI_v3'
bowtie_db_taxon_map_file = db_path + 'GM_BSI_v3_taxon_map.txt'
mothur_db_fa = db_path + 'mothur/silva_132_v3.fa'
mothur_db_tax = db_path + 'mothur/silva_132_v3.2.tax'

spades = bin_path + 'spades.py'
bowtie2 = bin_path + 'bowtie2'
cutadapt = bin_path + 'cutadapt'
trimmomatic = bin_path + 'trimmomatic'

cdhit = cmd_path + 'cd-hit'
mothur = cmd_path + 'mothur'
usearch = cmd_path + 'usearch11'

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
# 过滤拼接文库低质量
A1_trime = os.path.join(ana_dir, 'A1.fastq.gz')
submit_Trimmomatic_SE.submit_Trimmomatic_SE(A1, 
                                            A1_trime, 
                                            trimmomatic, 
                                            threads)
# 获得连接文库的R1和R2的UMI
logging.info(' UMI pairing start')
LR1=SeqIO.to_dict(SeqIO.parse(umi1_fq,'fastq'))
LR2=SeqIO.to_dict(SeqIO.parse(umi2_fq,'fastq'))
L_UMI_INFO_FILE = os.path.join(ana_dir,File_Tag+'.L.umi.INFO.txt')
L_UMI_ID_FILE = os.path.join(ana_dir,File_Tag+'.L.umiID.INFO.txt')
SEQ2UMI_FILE = os.path.join(ana_dir,File_Tag+'.seq2umi.txt')

df_L_UMI_INFO = get_L_UMI.get_L_UMI(LR1,LR2)
LR1 = {}
LR2 = {}
df_L_UMI_INFO.to_csv(L_UMI_INFO_FILE,sep='\t')

L_seq = df_L_UMI_INFO['counts_of_paire'].sum()
df_tmp = df_L_UMI_INFO[~(df_L_UMI_INFO['u1'].str.contains('\*'))]
L_R1_UMI = df_tmp['u1'].unique().shape[0]
L_R1_UMI_seq = df_tmp.drop_duplicates(['u1'])['u1_reads'].sum()
df_tmp = df_L_UMI_INFO[~(df_L_UMI_INFO['u2'].str.contains('\*'))]
L_R2_UMI = df_tmp['u2'].unique().shape[0]
L_R2_UMI_seq = df_tmp.drop_duplicates(['u2'])['u2_reads'].sum()

logging.info(' No. of L seqs:\t{}'.format(L_seq))
logging.info(' No. of L seqs found UMI(R1):\t{}'.format(L_R1_UMI_seq))
logging.info(' No. of UMIs in L data(R1):\t{}'.format(L_R1_UMI))
logging.info(' No. of L seqs found UMI(R2):\t{}'.format(L_R2_UMI_seq))
logging.info(' No. of UMIs in L data(R2):\t{}'.format(L_R2_UMI))
sta_list.append(['1 No. of L seqs:',L_seq])
sta_list.append(['2 No. of L seqs found UMI(R1):',L_R1_UMI_seq])
sta_list.append(['3 No. of UMIs in L data(R1):',L_R1_UMI])
sta_list.append(['4 No. of L seqs found UMI(R2):',L_R2_UMI_seq])
sta_list.append(['5 No. of UMIs in L data(R2):',L_R2_UMI])

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
    
tmp1 = df_L_UMI_INFO['umi_id'].unique().shape[0]
tmp2 = df_L_UMI_INFO.shape[0]
tmp3 = df_L_UMI_INFO['u1'].unique().shape[0]
tmp4 = df_L_UMI_INFO['u2'].unique().shape[0]
logging.info(' No. of UMI ID in L data:\t{}'.format(tmp1))
logging.info(' No. of UMI pairs in L data:\t{}'.format(tmp2))
logging.info(' No. of paired UMIs in L data(R1):\t{}'.format(tmp3))
logging.info(' No. of paired UMIs in L data(R2):\t{}'.format(tmp4))
sta_list.append(['6 No. of UMI ID in L data:',tmp1])
sta_list.append(['7 No. of UMI pairs in L data:',tmp2])
sta_list.append(['8 No. of paired UMIs in L data(R1):',tmp3])
sta_list.append(['9 No. of paired UMIs in L data(R2):',tmp4])
# 获取拼接文库中seq_id与umi对应关系
seq2umi,tmp_list = corresponding_seq2umi.corresponding_seq2umi(A2_cut_adapter_info,
                                                               A2_umi_file_info,
                                                               umi2ID)
sta_list.extend(tmp_list)
pd.DataFrame.from_dict(seq2umi,orient='index',columns=(['umi'])).to_csv(SEQ2UMI_FILE,sep='\t')
#分组
umi_seq_dir = os.path.join(ana_dir,'umi_seq')
mkdir.mkdir(umi_seq_dir)
logging.info(' write to umi paired seq file start')
(umi_counts,umi_id_counts,umi2seq) = group_umi_seq.queu_group_umi_seq(A1_trime,
                                                                      umi2ID,
                                                                      umi_seq_dir,
                                                                      seq2umi,
                                                                      threads
                                                                     )
logging.info(' write to umi paired seq file end')
#输出umi统计信息
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
logging.info(' No. of A reads after trimme:\t{}'.format(umi2seq[-1][-3]))
logging.info(' No. of UMI ID appeared in A data:\t{}'.format(umi2seq[-1][-1]))
sta_list.append(['11 No. of A reads after trimme:',umi2seq[-1][-3]])
sta_list.append(['17 No. of UMI ID appeared in A data:',umi2seq[-1][-1]])
#序列与umiID的关系图
umi2seq_fig = os.path.join(result_dir,File_Tag+'.umi2seq.pdf')
umi2seq_csv = os.path.join(result_dir,File_Tag+'.umi2seq.csv')
plot_seq2umi.plot_seq2umi(np.array(umi2seq),umi2seq_fig)
pd.DataFrame(np.array(umi2seq),
             columns=(['umi_seq','seq','umi','umi_id'])).to_csv(
                     umi2seq_csv,sep='\t')
#组装
#umi_id_fastq_file_list = os.path.join(result_dir,File_Tag+'.umi.fastq.list.txt')
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
    tmp_list.append(assemble_fastq_from_list.assemble(files,
                                                      spades_path_tmp,
                                                      spades_log_tmp,
                                                      spades,threads))
assmble_suc = [i[0] for i in tmp_list]
logging.info(' No.of Contigs assembled:\t{}'.format(sum(assmble_suc)))
sta_list.append(['18 No. of Contigs assembled:',sum(assmble_suc)])
if sum(assmble_suc) == 0:
    sys.exit('No Contigs!')
#组装后contig处理
#定义文件名称
merge_fa = os.path.join(ana_dir, File_Tag+'.merged.fasta')
# merge_fa_info = os.path.join(ana_dir, File_Tag+'.merged.fasta.info')
merge_trim_fa1 = os.path.join(ana_dir, File_Tag+'.merged.trim1.fasta')
merge_trim_log1 = os.path.join(ana_dir, File_Tag+'.merged.trim1.log')
merge_trim_fa2 = os.path.join(ana_dir, File_Tag+'.merged.trim2.fasta')
merge_trim_log2 = os.path.join(ana_dir, File_Tag+'.merged.trim2.log')
merge_filter_fa = os.path.join(ana_dir, File_Tag+'.merged.filter.fasta')
merge_filter_cluster_fa = os.path.join(ana_dir, File_Tag+'.merged.filter.rep.fa')
clust_fa = os.path.join(ana_dir, File_Tag+'.clust.fasta')
clust_fa_clstr_table = os.path.join(ana_dir, File_Tag+'.clust.clstr.info')
clust_uchimeout = os.path.join(ana_dir, File_Tag+'.clust.uchimeout.out.txt')
clust_ch_fa = os.path.join(ana_dir, File_Tag+'.clust.ch.fasta')
clust_nonch_fa = os.path.join(ana_dir, File_Tag+'.clust.nonch.fasta')
final_fa = os.path.join(result_dir, 'final.fasta')
final_tab = os.path.join(result_dir, 'final.size.tab.txt')
ID_info = os.path.join(result_dir, 'final.id.info')

#合并组装后的contigs
assemble_list = []
for tmpdir in os.listdir(spades_path):
    for d in os.listdir(os.path.join(spades_path, tmpdir)):
        contig_file = os.path.join(spades_path, tmpdir, d, 'contigs.fasta')
        if os.path.exists(contig_file):
            tmp = [d,contig_file]
            assemble_list.append(tmp)
df_merge_fa = merge_contigs.merge_contigs(assemble_list,merge_fa)
# df_merge_fa.to_csv(merge_fa_info,sep='\t')
#去除contig中引物序列
submit_cutadapt.submit_cutadapt(merge_fa, 
                                merge_trim_fa1, 
                                merge_trim_log1,
                                file_primer_rc,
                                'a',
                                cutadapt,
                                threads)
submit_cutadapt.submit_cutadapt(merge_trim_fa1,
                                merge_trim_fa2,
                                merge_trim_log2,
                                file_primer,
                                'g',
                                cutadapt,
                                threads)
#长度过滤
#大于minlength小于maxlength
tmp = cut_fa_by_len.cut_fa_by_len(merge_trim_fa2,
                                  merge_filter_fa,
                                  args.minlength,
                                  args.maxlength)
logging.info(' No. of Contigs after length filtered({}-{}):\t{}'.format(args.minlength, args.maxlength, tmp))
sta_list.append(['19 No. of Contigs after length filtered('+str(args.minlength)+'-'+str(args.maxlength)+'bp):',tmp])

#聚类
#100%相似度 cd-hit
submit_cdhit.submit_cdhit(merge_filter_fa, 
                          merge_filter_cluster_fa, 
                          1, 
                          cdhit, 
                          threads)
#去除同一umiID中未聚类在同一组的contig
consensus_seq_count, rep_seq_tab = get_consensus_seq_from_cdhit.get_consensus_seq_from_cdhit(
    merge_filter_cluster_fa+'.clstr', 
    merge_filter_fa, 
    clust_fa)
sta_list.append(['20 No. of Contigs after clust filter:',consensus_seq_count])
sta_list.append(['21 No. of Cluster:',len(rep_seq_tab)])
pd.DataFrame(rep_seq_tab,
             columns=(['clust_rep_id','size','seq_id'])).to_csv(clust_fa_clstr_table, 
                                                              sep='\t',
                                                              index=False)
#去除嵌合
remove_chimeras_uchime.remove_chimeras_uchime(clust_fa,
                                              clust_uchimeout,
                                              clust_ch_fa,
                                              clust_nonch_fa,
                                              usearch)
# 重新命名序列
dict_update_id = update_seqID.update_seqID(clust_nonch_fa, 
                                           final_fa, 
                                           File_Tag)
df_update_id = pd.DataFrame.from_dict(dict_update_id, 
                                      orient='index', 
                                      columns=(['contig_id']))
df_update_id['umi_id'] = pd.to_numeric(df_update_id['contig_id'].str.split('_').str[0])
df_ID_info = pd.merge(df_update_id.reset_index().rename(columns={'index':'ID'}), 
                      df_merge_fa.reset_index().rename(columns={'index':'umiID'}), 
                      left_on='umi_id', 
                      right_on='umi_id', 
                      how='left')
# 生成cluster size表格
tab_list = []
total_size = 0
total_cluster = 0
for i in SeqIO.parse(final_fa, 'fasta'):
    size = re.search('size=(\d+)', i.description).group(1)
    tab_list.append([i.id, size])
    total_size += int(size)
    total_cluster += 1
pd.DataFrame(tab_list, columns=(['id', File_Tag])).set_index('id').to_csv(final_tab, sep='\t')
sta_list.append(['22 No. of Cluster after remove Chimerias:',total_cluster])
sta_list.append(['23 No. of Contigs after remove Chimerias:',total_size])
# 生成统计表
sta_file = os.path.join(result_dir, 'sta.txt')
pd.DataFrame(sta_list,
             columns=(['description', File_Tag])).to_csv(sta_file, 
                                                         sep='\t', 
                                                         index=False)

#物种分类
submit_mothur.submit_mothur(final_fa, mothur_db_fa, mothur_db_tax, mothur, threads)
tax_file = (re.search(r'(.*)\.fasta',final_fa).group(1) + 
              re.search(r'(\..*)\.tax',mothur_db_tax).group(1) + 
              '.wang.taxonomy')
tax_ratio_file = os.path.join(result_dir, 'tax.ratio.xlsx')
tax_reads_file = os.path.join(result_dir, 'tax.reads.xlsx')
get_tax_from_mothur_with_level.get_tax_from_mothur_with_level(tax_file,
                                                              final_tab,
                                                              tax_ratio_file,
                                                              tax_reads_file)
if args.group:
    target_dir = os.path.join(result_dir, File_Tag+'_umiID_Reads')
    mkdir.mkdir(target_dir)
    file_path = {}
    for tmpdir in os.listdir(umi_seq_dir):
        path = os.path.join(umi_seq_dir, tmpdir)
        for file in os.listdir(path):
            f = os.path.splitext(file)[0]
            file_path[f] = os.path.join(path, file)
    for u in df_ID_info['umiID'].to_list():
        t_dir = os.path.join(target_dir, u[0:2])
        mkdir.mkdir(t_dir)
        os.system(' '.join(['/usr/bin/mv', file_path[u], t_dir]))
    os.system(' '.join(['/usr/bin/tar', 'czfP', target_dir+'.tar.gz', target_dir]))
    os.system(' '.join(['/usr/bin/rm','-r',  target_dir]))
    df_ID_info.to_csv(ID_info, sep='\t')
    
if not args.debug:
    shutil.rmtree(ana_dir)
    df_ID_info.to_csv(ID_info, sep='\t')
logging.info(' All Finishied!')


