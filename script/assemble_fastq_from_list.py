#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:56:20 2019
序列组装
@author: yk
"""

import argparse
import logging
import multiprocessing
import os
import re
from script import mkdir


def spades_sub(file_list):
    '''提交spades组装任务
    参数：
        list: [file, outdir, logfile, spades_path]
            file: 用于组装的fastq序列文件
            outdir: 组装结果输出目录
            logfile: 组装日志输出文件
            spades_path: spades.py程序路径
    返回：
        程序执行的返回值（0,1）
    '''
    file = file_list[0]
    outdir = file_list[1]
    logfile = file_list[2]
    spades_path = file_list[3]
    cmd = spades_path + ' --s1 ' + file + ' -o ' + outdir + ' -t 1 ' \
          + ' 1>' + logfile + ' 2>' + logfile
    tmp_sys = os.system(cmd)
    rm_cmd = ('/bin/rm -rf outdir/assembly_graph* '
              'outdir/before_rr.fasta '
              'outdir/contigs.paths '
              'outdir/corrected '
              'outdir/dataset.info '
              'outdir/input_dataset.yaml '
              'outdir/misc '
              'outdir/K* '
              'outdir/params.txt '
              'outdir/tmp '
              'outdir/scaffolds.* '
              'outdir/spades.log '
              'outdir/warnings.log '
              'outdir/configs '
              )
    try:
        tmp_sys_rm = os.system(re.sub('outdir', outdir, rm_cmd))
    except OSError as e:
        logging.info('rm error {} - {}'.format(e.filename, e.strerror))
    return tmp_sys


def assemble(file_list, out_dir, log_dir, spades_path, threads=1, ):
    '''多进程提交序列组装任务
    参数：
        file_list: 用于组装的序列文件列表
        out_dir: 组装结果输出目录
        log_dir: 组装日志输出目录
        spades_path: spades.py程序路径
        threads: 进程数
    返回：
        ass_success: 组装成功(结果中有contigs.fasta)数量
        ass_failed: 组装失败数量
    '''
    # tmp_file_list=list(file_list)[:]
    ass_success = 0
    ass_failed = 0
    mkdir.mkdir(out_dir)
    mkdir.mkdir(log_dir)
    queu_file = []
    for f in file_list:
        f_name = os.path.splitext(os.path.basename(f))[0]
        spades_tmp_dir = os.path.join(out_dir, f_name)
        logfile = os.path.join(log_dir, f_name)
        queu_file.append([f, spades_tmp_dir, logfile, spades_path])
    pool = multiprocessing.Pool(processes=threads)
    spades_return = pool.map_async(spades_sub, queu_file)
    spades_return.wait()

    if spades_return.ready():
        for get in spades_return.get():
            if get == 0:
                ass_success += 1
            else:
                ass_failed += 1

    pool.close()
    pool.join()
    return ass_success, ass_failed


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',
                        level='INFO')
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="assemble fastq multiproced using spades.py from fastq dir")
    parse.add_argument('-i', '--indir',
                       required=True,
                       help='fastq dir which has nothing but fastq file in it')
    parse.add_argument('-o', '--outdir',
                       required=True,
                       help='out dir of assemble result')
    parse.add_argument('-l', '--logdir',
                       default='',
                       help='out log dir of assemble log.default=outdir_log')
    parse.add_argument('-t', '--threads',
                       default=4,
                       type=int,
                       help='threads for assemble fastq')
    parse.add_argument('-p', '--spadespath',
                       default='/Bioinfo/bin/spades.py',
                       help='spades.py absolute path. ')
    args = parse.parse_args()

    if args.logdir:
        logdir = args.logdir
    else:
        logdir = args.outdir + '_log'

    files = map(lambda x: args.indir + '/' + x, os.listdir(args.indir))
    logging.info(' assemble {} files using {} by {} threads:\t' \
                 .format(len(os.listdir(args.indir)), args.spadespath, args.threads))
    ass_success, ass_failed = assemble(files,
                                       args.outdir,
                                       logdir,
                                       args.spadespath,
                                       args.threads)
    logging.info(' assemble sucessed {}:\t'.format(ass_success))
    logging.info(' assemble failed {}:\t'.format(ass_failed))
