#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: yk
"""

import argparse
import logging
import os
import subprocess


def submit_mothur(infa, taxdb, tax, mothur, threads):
    '''物种分类
    参数：
        infa: 输入fasta文件
        taxdb: mothur数据库
        tax: mothur数据库分类文件
        mothur: mothur路径
        threads: 进程数
    返回：
        无
    '''
    cmd = mothur
    cmd = cmd + ' "#classify.seqs(fasta='
    cmd = cmd + infa + ',template='
    cmd = cmd + taxdb + ',taxonomy='
    cmd = cmd + tax + ',processors='
    cmd = cmd + str(threads) + ')"'
    logging.info(' %s' % cmd)
    try:
        status = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        print(status)
    except Exception as e:
        logging.info(' %s' % cmd)


if __name__ == '__main__':
    db_dir = '/Bioinfo/Database/Silva_132_v3/mothur/'
    parser = argparse.ArgumentParser(
        description='this is for taxonomy of fasta using mothur classify.seqs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--in_fasta',
                        help='input fasta file',
                        required=True)
    parser.add_argument('-db', '--taxdb',
                        help='tax database fasta file',
                        default=os.path.join(db_dir, 'silva_132_v3.fa'))
    parser.add_argument('-tax', '--tax',
                        help='tax database tax file',
                        default=os.path.join(db_dir, 'silva_132_v3.2.tax'))
    parser.add_argument('-t', '--threads',
                        help='threads for bowtie.',
                        default=1, type=int)
    parser.add_argument('-path', '--mothur_path',
                        help='mothur absolute path',
                        default='/Bioinfo/bin/mothur')
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)
    submit_mothur(args.in_fasta, args.taxdb, args.tax, args.mothur_path, args.threads)
