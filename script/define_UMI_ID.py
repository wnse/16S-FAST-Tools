#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:48:08 2019
根据连接文库序列的umi，定义umiID。
1.排序，确定数量最多的配对关系；
2.在剩余的配对关系中删除：
    已经确定的umiID的umi2对应的所有其他umi配对中，该umi2为最高数量的umi1，除非该配对关系数量大于第一配对关系数量的特定倍数；
    已经确定的umiID的umi1对应的所有其他umi配对中，该umi1为最高数量的umi2，除非该配对关系数量大于第一配对关系数量的特定倍数；
3.如果配对关系数量小于指定的最小值，配对结束；
4.如果已无配对，配对结束；
@author: yk
"""
import argparse
import logging
import pandas as pd


def define_UMI_ID(UMI_dict, UMI1_dict, UMI2_dict, cutoff=0.5, min_counts=1):
    ''' 根据连接文库序列的umi，定义umiID
    参数：
        UMI_dict: 连接文库所有umi配对关系及数量（type=dict）
        UMI1_dict: 连接文库所有umi1对应的umi2及数量
        UMI2_dict: 连接文库所有umi2对应的umi1及数量
        cutoff: 可以作为umiID的配对关系所需要大于第一配对关系数量的最小比例
        min_counts: 配对关系的最小数量
    返回：
        out_list: [umiID, umi1, umi2, counts]
    '''
    out_list = []
    UMI_list = sorted(UMI_dict.items(), key=lambda UMI_dict: UMI_dict[1], reverse=True)
    for idx_count in UMI_list:
        if len(UMI1_dict) == 0:
            break
        (u1, u2) = idx_count[0].split('_')
        if u1 in UMI1_dict.keys():
            if len(UMI1_dict[u1]) == 0:
                del UMI1_dict[u1]
                continue
            if u2 not in UMI2_dict.keys():
                del UMI1_dict[u1][u2]
                continue
            if UMI1_dict[u1][u2] < min_counts:
                logging.info(' exit for {} counts {} < {}' \
                             .format(idx_count[0], UMI1_dict[u1][u2], min_counts))
                break
            out_list.append([idx_count[0], u1, u2, UMI1_dict[u1][u2]])
            list1 = list(UMI2_dict[u2].keys())
            list2 = list(UMI1_dict[u1].keys())
            max_counts = UMI1_dict[u1][u2]
            del UMI1_dict[u1]
            del UMI2_dict[u2]
            for deu in list1:
                if deu in UMI1_dict.keys():
                    if max(UMI1_dict[deu], key=UMI1_dict[deu].get) == u2:
                        if UMI1_dict[deu][u2] >= cutoff * max_counts:
                            out_list.append([idx_count[0], deu, u2, UMI1_dict[deu][u2]])
                        del UMI1_dict[deu]
            for deu in list2:
                if deu in UMI2_dict.keys():
                    if max(UMI2_dict[deu], key=UMI2_dict[deu].get) == u1:
                        if UMI2_dict[deu][u1] >= cutoff * max_counts:
                            out_list.append([idx_count[0], u1, deu, UMI2_dict[deu][u1]])
                        del UMI2_dict[deu]
    return out_list


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level='INFO')
    parse = argparse.ArgumentParser(description='get max pair UMI (UMI_ID) from L UMI INFO file',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('L_UMI_INFO_FILE',
                       help='L UMI INFO file(csv) from get_L_UMI.py')
    parse.add_argument('-c', '--cutoff',
                       help='threshold for another umi paire in a UMI_ID,',
                       default=0.5,
                       choices=([i / 10 for i in range(1, 10)]))
    parse.add_argument('-m', '--min_count',
                       help='min paires for umis to be an UMI_ID,',
                       default=1, type=int)
    args = parse.parse_args()

    df = pd.read_csv(args.L_UMI_INFO_FILE, sep='\t', index_col=0)
    df = df[~ ((df['u1'].str.contains('\*')) | (df['u2'].str.contains('\*')))]

    umi = {}
    umi1 = {}
    umi2 = {}
    for i in list(df.index):
        u1 = str(df.loc[i, 'u1'])
        u2 = str(df.loc[i, 'u2'])
        u = u1 + '_' + u2
        c = int(df.loc[i, 'counts_of_paire'])
        umi[u] = c
        umi1[u1] = umi1.get(u1, {})
        umi1[u1][u2] = c
        umi2[u2] = umi2.get(u2, {})
        umi2[u2][u1] = c

    tmp_list = define_UMI_ID(umi, umi1, umi2, args.cutoff, args.min_count)
    df = pd.DataFrame(tmp_list)
    df.columns = (['umi_id', 'u1', 'u2', 'counts_of_paire'])
    df.to_csv('L_UMI_ID.txt', sep='\t', index=False)
