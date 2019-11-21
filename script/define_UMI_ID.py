#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:48:08 2019

@author: yk
"""
import argparse
import logging
import pandas as pd

def define_UMI_ID(UMI_dict,UMI1_dict,UMI2_dict,cutoff=0.5,min_counts=1):
    ''' 
    1.排序UMI_dict中counts最多的idx；
    2.在剩余的idx中删除：
        已经确定的的UMI对的UMI2对应的所有其他UMI对中，该UMI2为最高数量的UMI1,除非该配对counts大于最多counts的cut_off倍数；
        已经确定的的UMI对的UMI1对应的所有其他UMI对中，该UMI1为最高数量的UMI2,除非该配对counts大于最多counts的cut_off倍数；
    3.如果配对counts小于min_counts，结束；
    4.如果无配对，循环结束；
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

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
    parse = argparse.ArgumentParser(description='get max pair UMI (UMI_ID) from L UMI INFO file')
    parse.add_argument('L_UMI_INFO_FILE',help='L UMI INFO file(csv) from get_L_UMI.py')
    parse.add_argument('-c','--cutoff',help='threshold for another umi paire in a UMI_ID, \
                        default=0.5,e.g.,more than 0.5*(first umi paire)',\
                        default=0.5,choices=([i/10 for i in range(1,10)]))
    parse.add_argument('-m','--min_count',help='min paires for umis to be an UMI_ID, \
                       default=1',default=1,type=int)
    args = parse.parse_args()
    
    df = pd.read_csv(args.L_UMI_INFO_FILE,sep='\t',index_col=0)
    df = df[~ ((df['u1'].str.contains('\*')) | (df['u2'].str.contains('\*')))]
    
    umi = {}
    umi1 = {}
    umi2 = {}
    for i in list(df.index):
        u1 = str(df.loc[i,'u1'])
        u2 = str(df.loc[i,'u2'])
        u = u1 + '_' + u2
        c = int(df.loc[i,'counts_of_paire'])
        umi[u] = c
        umi1[u1] = umi1.get(u1,{})
        umi1[u1][u2] = c
        umi2[u2] = umi2.get(u2,{})
        umi2[u2][u1] = c
    
    
    tmp_list = define_UMI_ID(umi,umi1,umi2,args.cutoff,args.min_count)
    df = pd.DataFrame(tmp_list)
    df.columns=(['umi_id','u1','u2','counts_of_paire'])
    df.to_csv('L_UMI_ID.txt',sep='\t',index=False)
    