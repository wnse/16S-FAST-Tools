import argparse
import os
import logging


def submit_cdhit(infa, outfa, identity, path, threads=1):
    '''序列聚类
    参数：
        infa: 输入fasta序列
        outfa: 输出fasta代表序列
        identity: 相似度
        path: cd-hit路径
        threads: 线程数
    返回：
        无
    '''
    cdhit = path
    str_join = ' '
    cmd = str_join.join([cdhit,
                         '-i', infa,
                         '-o', outfa,
                         '-c', str(identity),
                         '-T', str(threads)])
    logging.info(' %s' % cmd)
    status = os.system(cmd)
    if status == 0:
        logging.info(' done')
    else:
        logging.info(' exit')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='this is for cluster of seq using cd-hit',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infasta',
                        help='input fasta',
                        required=True)
    parser.add_argument('-o', '--outfasta',
                        help='output fastq',
                        required=True)
    parser.add_argument('-t', '--threads',
                        help='threads for computing.',
                        default=1,
                        type=int)
    parser.add_argument('-c', '--identity',
                        help='sequence identity threshold',
                        default=1,
                        type=float)
    parser.add_argument('-path', '--path',
                        help='cd-hit absolute path',
                        default='/root/anaconda3/bin/cd-hit')
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)
    submit_cdhit(args.infasta, args.outfasta, args.identity, args.path, args.threads)
