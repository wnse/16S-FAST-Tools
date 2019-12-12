import argparse
import logging
import os

def submit_Trimmomatic_SE(in_file,out_file,path,threads=1):
    str_join=' '
    trimmomatic = path
    cmd=str_join.join([trimmomatic,
                       'SE -phred33 ',
                       '-threads',str(threads),
                       in_file,
                       out_file,
                       'SLIDINGWINDOW:10:30',
                       'MINLEN:100',])
#                        'ILLUMINACLIP:'+adapter+'2:30:10'])
    logging.info(' %s'%cmd) 
    status=os.system(cmd)
    if status==0:
        logging.info(' done')
    else:
        logging.info(' exit')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fq','--fastq',
                        help='input fastq',required=True)
    parser.add_argument('-o','--out_fastq',
                        help='output fastq',required=True)
    parser.add_argument('-t','--threads',
                        help='threads for computing.',
                        default=1,type=int)
    parser.add_argument('-path','--path',
                        help='cutadapt absolute path,',
                        default='/root/anaconda3/bin/trimmomatic')
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
    submit_Trimmomatic_SE(args.fastq,
                          args.out_fastq,
                          args.path,
                          args.threads,
                         )
