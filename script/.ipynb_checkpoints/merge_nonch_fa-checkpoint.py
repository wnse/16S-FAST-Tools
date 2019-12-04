#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:27:39 2019

@author: yk
"""
import argparse
import re
import pandas as pd
from Bio import SeqIO

def merge_nonch_fa(ap):
# merge nonch.fa
    fa_list = [i for i in ap.fasta.split(':') if(len(str(i))!=0)]
    outfa = ap.out_fasta
    asvtab = ap.out_tabble
    tag = ap.tag
    
    seq_size = {}
    for fa in falist:
        if os.path.exists(fa):
            print(fa)
        else:
            print(fa + ' does not exists')
            continue
        seqs = SeqIO.parse(fa,'fasta')
        for seq in seqs:
            key = str(seq.seq)
            name = re.search(r'^(.*)_\d+$',str(seq.id)).group(1)
            size = re.search(r'\s+(.*);size=(\d+)',str(seq.description)).group(2)
            #name = match.group(1)
            #size = match.group(2)
            seq_size[key] = seq_size.get(key,{})
            seq_size[key][name] = seq_size[key].get(name,0) + int(size)
            
    df = pd.DataFrame.from_dict(seq_size).T
    df['sum'] = df.sum(axis=1)
    df.sort_values(by='sum',inplace=True,ascending=False)
    print(df.describe())
    df.reset_index(inplace=True)
    out = open(outfa,'w')
    for i in df.index:
        print('>' + tag + '_' + str(i) + ' asv;size=' + str(df.loc[i,'sum']),file = out)
        print(str(df.loc[i,'index']),file = out)
    out.close()
    df.drop(['index','sum'],axis=1,inplace=True)
    df.fillna(0).to_csv(asvtab,sep='\t')


def asvfa2df(fa):
    seqs = SeqIO.parse(fa,'fasta')
    seq_dict={}
    for seq in seqs:
        s = str(seq.seq)
        i = str(seq.id)
        name = re.search(r'^(.*)_\d+$',str(seq.id)).group(1)
        size = re.search(r'\s+(.*);size=(\d+)',str(seq.description)).group(2)
        seq_dict[s] = seq_dict.get(s,{})
        seq_dict[s][name] = size
        seq_dict[s][name+'_id'] = i
    df = df = pd.DataFrame.from_dict(seq_dict).T
    return df,name

def merge_nonch_fa_sample2batch(ap):
    sample_fa = ap.sample_fa
    batch_fa = ap.batch_fa
    batch_asv_tab = ap.batch_tab
    outfa = ap.out_fa
    outtab = ap.out_tab
    sample_tax_mothur = ap.sample_tax
    batch_tax_mothur = ap.batch_tax
    outtax = ap.out_tax
    tag = ap.tag

    df_batch_seq,batch_name = asvfa2df(batch_fa)
    df_sample_seq,sample_name = asvfa2df(sample_fa)
    df = pd.merge(df_batch_seq,df_sample_seq,left_index=True,right_index=True,how='outer') #合并序列
    df['sum'] = df[[batch_name,sample_name]].apply(pd.to_numeric).sum(axis=1) #计算合并后的size
    df = df.sort_values(by='sum',ascending=False).reset_index() #根据size排序并重新index
    df['tmp_id'] = df.index
    df['ID']='asv_'
    df['ID']=df['ID'].str.cat(df['tmp_id'].apply(str)) #形成新的ID名称
    df_write = df[['ID','sum','index']]
    write_df_to_fa(df_write,outfa)
    
    df_asv = pd.read_csv(batch_asv_tab,index_col=0,sep='\t')
    df_asv['asv_id']='asv_'
    df_asv['asv_id']=df_asv['asv_id'].str.cat(pd.Series(list(df_asv.index)).apply(str)) #加上 asv_ 字样，与fa的ID一致
    df_asv = pd.merge(df_asv,df[['tmp_id',batch_name+'_id',sample_name]],\
                      left_on='asv_id',right_on=batch_name+'_id',how='outer') #利用原有ID进行asv tab的合并
    df_asv = df_asv.sort_values(by='tmp_id').reset_index() #根据新的ID进行排序
    df_asv['new_id']='asv_'
    df_asv['new_id']=df_asv['new_id'].str.cat(pd.Series(list(df_asv.index)).apply(str)) #形成新的ID名称
    df_asv.drop(['asv_id','tmp_id','index'],axis=1,inplace=True)
    df_asv.set_index('new_id',inplace=True)
    df_asv.fillna(0).to_csv(outtab,sep='\t')
    
    if batch_tax_mothur and sample_tax_mothur:
        df_batch_tax = pd.read_csv(batch_tax_mothur,sep='\t',header=None,index_col=0)
        df_batch_tax.columns=(['batch_tax'])
        df_batch_tax['batch_tax_n'] = rm_score_from_tax(df_batch_tax['batch_tax'])
        
        df_sample_tax = pd.read_csv(sample_tax_mothur,sep='\t',header=None,index_col=0)
        df_sample_tax.columns=(['sample_tax'])
        df_sample_tax['sample_tax_n'] = rm_score_from_tax(df_sample_tax['sample_tax'])
        
        df_merge_tax = pd.merge(df[['ID',batch_name+'_id',sample_name+'_id']].fillna('-'),df_batch_tax,\
                                left_on=batch_name+'_id',right_index=True,how='left')
        df_merge_tax = pd.merge(df_merge_tax.fillna('-'),df_sample_tax,\
                                left_on=sample_name+'_id',right_index=True,how='left')
        df_merge_tax.fillna('-',inplace=True)
        df_merge_tax['consist'] = list(map(lambda x,y:x==y,list(df_merge_tax['batch_tax_n']),\
                                           list(df_merge_tax['sample_tax_n'])))

#        df_unconsist = df_merge_tax[~((df_merge_tax['batch_tax']=='-')\
#                                      |(df_merge_tax['sample_tax']=='-'))\
#                                    & df_merge_tax['consist'] == 'False']
#        if not df_unconsist.empty:
#            print (df_unconsist.drop(['batch_tax_n','sample_tax_n'],axis=1))
            
        keep = 'batch_tax'
        other = 'sample_tax'
        keep_n = df_merge_tax[df_merge_tax[keep]=='-']
        df_merge_tax.loc[keep_n.index,[keep]] = keep_n[other]
        df_merge_tax[['ID',keep]].to_csv(outtax,sep='\t',header=None,index=False)

def rm_score_from_tax(df):
    return df.str.replace('\(\d+\)','')

def write_df_to_fa(df_input,outfa):
    out = open(outfa,'w')
    for i in df_input.index:
        print('>' + str(df_input.loc[i][0]) + ' asv;size=' + str(df_input.loc[i][1]),file = out)
        print(str(df_input.loc[i][2]),file = out)
    out.close()

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.add_subparsers(help='merge sample nonch fasta to a fasta and table or \
                                 merge sample to batch',)#dest='subparser_name')
    
    args1 = args.add_parser('samples',help='merge each sample nonch fasta to a batch fasta and table')
    args1.add_argument('-fa','--fasta',help='fasta files joined by ":", eg "test1.fa:test2.fa:..."',required=True,type=str)
    args1.add_argument('-outfa','--out_fasta',help='output merged fasta file.',default='total_merged.fa',type=str)
    args1.add_argument('-table','--out_table',help='output merged asv table.',default='total_merged_tab.txt',type=str)   
    args1.add_argument('-tag','--tag',help='fasta sequences description tag',default='asv')   
    args1.set_defaults(func=merge_nonch_fa)
    
    args2 = args.add_parser('sample2batch',help='merge a samle nonch fasta to a batch fasta and table.')
    args2.add_argument('-sfa','--sample_fa',help='sample nonch fasta file.',required=True,type=str)
    args2.add_argument('-stax','--sample_tax',help='sample mothur tax file.')
    args2.add_argument('-bfa','--batch_fa',help='batch fasta.',required=True,type=str)
    args2.add_argument('-btab','--batch_tab',help='batch table.',required=True,type=str)
    args2.add_argument('-btax','--batch_tax',help='batch mothur tax file.')
    args2.add_argument('-ofa','--out_fa',help='output merged fasta.',default='total_merged.fa')
    args2.add_argument('-otab','--out_tab',help='output merged table.',default='total_merged_tab.txt')
    args2.add_argument('-otax','--out_tax',help='output merged tax file',default='total_merged_tax.txt')
    args2.add_argument('-tag','--tag',help='fasta sequences description tag',default='asv')
    args2.set_defaults(func=merge_nonch_fa_sample2batch)
    
    ap = parser.parse_args()
    print(ap)
    ap.func(ap)

#    fa_list = [i for i in args.fasta.split(':') if(len(str(i))!=0)]
#    merge_nonch_fa(fa_list,args.out_fasta,args.out_table
