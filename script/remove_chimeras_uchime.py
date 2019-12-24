import os
import sys
import argparse


def remove_chimeras_uchime(fa, uchimeout, chf, nonchf, usearch_path):
    '''去除嵌合序列
    参数：
        fa: 输入fasta文件
        uchimeout: 去除嵌合过程文件
        chf: 嵌合序列
        nonchf: 非嵌合序列
        usearch_path: usearch路径
    返回：
        无
    '''
    cmd = ' '.join([usearch_path,
                    '-uchime3_denovo', fa,
                    '-uchimeout', uchimeout,
                    '-chimeras', chf,
                    '-nonchimeras', nonchf])
    os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-fa', '--fasta',
                        help='input fasta file',
                        required=True)
    parser.add_argument('-out', '--out',
                        help='uchime out txt. default=input_fasta.uchimeout.out.txt',
                        default='')
    parser.add_argument('-chfa', '--chimeras_fasta',
                        help='output chimeras fasta file. default=input_fasta.ch.fa',
                        default='')
    parser.add_argument('-nonchfa', '--nonchimeras_fasta',
                        help='output non-chimeras fasta file. default=input_fasta.nonch.fa',
                        default='')
    parser.add_argument('-usearch', '--usearch_path',
                        help='absolute path of usearch software. default=/Bioinfo/usearch11',
                        default='/Bioinfo/bin/usearch11')
    args = parser.parse_args()

    out = ''
    ch = ''
    nonch = ''
    usearch = ''
    if args.out:
        out = args.out
    else:
        out = args.fasta + '.uchimeout.out.txt'

    if args.chimeras_fasta:
        ch = args.chimeras_fasta
    else:
        ch = args.fasta + '.ch.fa'

    if args.nonchimeras_fasta:
        nonch = args.nonchimeras_fasta
    else:
        nonch = args.fasta + '.nonch.fa'

    if os.path.exists(args.usearch_path):
        usearch = args.usearch_path
    else:
        sys.exit(args.usearch_path + ' does not exists')
    remove_chimeras_uchime(args.fasta, out, ch, nonch, args.usearch_path)
