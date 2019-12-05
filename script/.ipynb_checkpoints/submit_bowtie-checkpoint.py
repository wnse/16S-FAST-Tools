#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 11:46:06 2019

@author: yk
"""

import argparse
import logging
import os

def submit_bowtie(in_file,out_file,log_file,threads,bowtie_path):
    '''
    bowtie alignment
    '''
    str_join=' '
    bowtie2 = path
    cmd=str_join.join([bowtie2,\
                       '-p', str(threads),\
                       '-x',bowtie_db,\
                       '-f',in_file,\
                       '-S',out_file,\
                       '1>',log_file,\
                       '2>',log_file])
    logging.info(' %s'%cmd)
    status=os.system(cmd)
    if status==0:
        logging.info(' done')
    else:
        logging.info(' exit')
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
                                     'this is for alignment of fasta \
                                     using bowtie2')
    parser.add_argument('-i','--in_fasta',help='input fasta file',required=True)
    parser.add_argument('-o','--out_sam',help='output sam file',required=True)
    parser.add_argument('-l','--log',help='output log file. default=out.log',\
                        default='')
    parser.add_argument('-t','--threads',help='threads for bowtie. default=1',default=1,type=int)
    parser.add_argument('-path','--bowtie2_path',help='bowtie2 absolute path. \
                        default=/Users/yk/anaconda3/envs/bioconda/bin/bowtie2',\
                        default='/Users/yk/anaconda3/envs/bioconda/bin/bowtie2')
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
    log_check = args.log if args.log else args.out + '.log'
    submit_bowtie(args.in_fasta,args.out_sam,log_check,args.threads,args.bowtie_path)    