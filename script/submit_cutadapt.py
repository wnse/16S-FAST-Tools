#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:24:34 2019

@author: yk
"""

import argparse
import logging
import os

def submit_cutadapt(in_file,out_file,log_file,adapter,g_a,path,threads=1,info=0):
    '''
    cutadapt remove fastq primer
    '''
    str_join=' '
    info_file = out_file + '.cutadapt.info.file'
    cutadapt = path
    if info:
        cmd=str_join.join([cutadapt,'-%s'%g_a,adapter,'-o',out_file,in_file,'--info-file',info_file,'>',log_file])
    else:
        cmd=str_join.join([cutadapt,'-%s'%g_a,adapter,'-o',out_file,in_file,'-j',str(threads),'>',log_file]) 
    logging.info(' %s'%cmd) 
    status=os.system(cmd)
    if status==0:
        logging.info(' done')
    else:
        logging.info(' exit')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=\
                                     'this is for trimming adapter of fasta \
                                     using cutadapt')
    parser.add_argument('-fq','--fastq',\
                        help='input fastq',required=True)
    parser.add_argument('-o','--out_fastq',\
                        help='output fastq',required=True)
    parser.add_argument('-l','--log_file',\
                        help='output log file. default=out_fastq.log',default='')
    parser.add_argument('-a','--adapter',\
                        help='adapter sequences',required=True)
    parser.add_argument('-d','--direction',\
                        help="g(5') or a(3') for trimming. default=a.",\
                        choices=['g','a'],default='a')
    parser.add_argument('-t','--threads',\
                        help='threads for computing. default=out_fastq.log',\
                        default=1,type=int)
    parser.add_argument('-i','--info',\
                        help='output info file named out_fastq.cutadapt.info.file,\
                        if set -i,threads not useful',\
                        action="store_true")
    parser.add_argument('-path','--path',\
                        help='cutadapt absolute path,\
                        default=/Bio/bin/cutadapt',\
                        default='/Bio/bin/cutadapt')
    
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
    log_file = args.log_file if args.log_file else args.out_fastq + '.log'
    if args.info:
        info_check = 1
    else:
        info_check = 0
    submit_cutadapt(args.fastq,\
                    args.out_fastq,\
                    log_file,\
                    args.adapter,\
                    args.direction,\
                    args.path,\
                    args.threads,
                    info_check)
