import argparse
from Bio import SeqIO
import pandas as pd


def update_seqID(infa, outfa, tag):
    '''重新命名序列
    参数：
        infa: 输入fasta文件
        outfa: 输出fasta文件
        tag: 序列标签
    返回：
        id_track: 序列ID对应关系（type=dict）
    '''
    n = 0
    id_track = {}
    seqs = []
    for rec in SeqIO.parse(infa, 'fasta'):
        new_id = tag + '_' + str(n)
        id_track[new_id] = rec.id
        rec.id = new_id
        rec.name = new_id
        rec.description = str(rec.description).split()[-1]
        n += 1
        seqs.append(rec)
    count = SeqIO.write(seqs, outfa, 'fasta')
    return id_track


if __name__ == '__main__':
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', '--infa', required=True)
    parse.add_argument('-o', '--outfa', required=True)
    parse.add_argument('-t', '--tag', required=True)
    args = parse.parse_args()
    tmp_dict = update_seqID(args.infa, args.outfa, args.tag)
    print(pd.DataFrame.from_dict(tmp_dict, orient='index', columns=(['original_id'])))
