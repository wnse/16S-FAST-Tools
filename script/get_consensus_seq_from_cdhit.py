import argparse
import logging
import re
import pandas as pd
from Bio import SeqIO


def get_consensus_seq_from_cdhit(clstr, fa, outfa):
    '''从cd-hit聚类结果中过滤一致性umiID
    参数：
        clstr: cd-hit聚类结果文件
        fa: 用于聚类的所有序列fa文件
        outfa: 输出结果文件
    返回：
        total_size: 过滤后的序列数量
        rep_seq_tab: 代表序列及对应clust的序列ID（type=pandas）
    '''
    # 读取每个id的所有序列对应的cluster
    seq2cl = {}
    for rec in open(clstr, 'rt'):
        if re.match(r'^>', rec):
            tmp = re.search(r'^>(.+)', rec).group(1)
        else:
            seq = re.search(r'>(.+)\.\.\.', rec).group(1)
            name, name_i = re.search('(.+)_(\d+)', seq).group(1, 2)
            seq2cl[name] = seq2cl.get(name, [])
            seq2cl[name].append([name_i, tmp])
    # 读取每条序列的长度
    seqLen = {}
    for rec in SeqIO.parse(fa, 'fasta'):
        seqLen[rec.id] = len(rec.seq)
    # 若id的所有序列分布在不同的cluster，舍去id
    # 反之，将选择该id的最长的一条序列加入cluster
    cl2seq = {}
    for name, cl in seq2cl.items():
        cl_s = ([i[1] for i in cl])
        if len(set(cl_s)) == 1:
            tmpLen = 0;
            tmpcl = ''
            for i in cl:
                seqname = name + '_' + i[0]
                if seqLen[seqname] > tmpLen:
                    tmpcl = seqname
                    tmpLen = seqLen[seqname]
            cl2seq[cl_s[0]] = cl2seq.get(cl_s[0], [])
            cl2seq[cl_s[0]].append(tmpcl)
        else:
            logging.info('unconsensus seq {} in cluster {}'.format(name, set(cl_s)))
    # 在每个cluster中选择最长序列为代表序列
    rep_seq_size = {}
    rep_seq_seqs = {}
    total_size = 0
    for cl, seqs in cl2seq.items():
        tmp = [seqLen[i] for i in seqs]
        rep_seq = seqs[tmp.index(max(tmp))]
        rep_seq_size[rep_seq] = len(seqs)
        total_size += len(seqs)
        rep_seq_seqs[rep_seq] = seqs
    # 按照cluster大小顺序输出代表序列
    rep_seq = []
    rep_seq_tab = []
    tmp_seq = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    n = 0
    for i, size in sorted(rep_seq_size.items(), key=lambda d: d[1], reverse=True):
        seq = tmp_seq[i]
        #         new_id = tag + '_' + str(n)
        #         seq.id = new_id
        seq.description = seq.description + ';size=' + str(size)
        rep_seq.append(seq)
        rep_seq_tab.append([seq.id, size, '|'.join(rep_seq_seqs[i])])
        n += 1
    count = SeqIO.write(rep_seq, outfa, 'fasta')
    return total_size, rep_seq_tab


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inclstr',
                        help='cd-hit clstr file',
                        required=True)
    parser.add_argument('-f', '--infasta',
                        help='fasta file used to cluster by cd-hit',
                        required=True)
    parser.add_argument('-o', '--outfasta',
                        help='output fasta',
                        required=True)
    parser.add_argument('-tab', '--outtab',
                        help='output cluster table. default= *.clstr.table',
                        default='')
    parser.add_argument('-name', '--name',
                        help='sample name',
                        default='cluster')
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',
                        level='INFO')
    sample_name = 'test'
    outtab = ''
    if args.outtab:
        outtab = args.outtab
    else:
        outtab = args.outfasta + '.clstr.table'
    total_size, rep_seq_tab = get_consensus_seq_from_cdhit(args.inclstr,
                                               args.infasta,
                                               args.outfasta)
    pd.DataFrame(rep_seq_tab,
                 columns=(['asv_id', 'size', 'seq_id'])).to_csv(outtab,
                                                                sep='\t',
                                                                index=False)
    logging.info(' Total Contig after filterd:\t{}'.format(total_size))
