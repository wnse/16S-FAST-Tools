import json
import os
import pandas as pd
import argparse
import time
import logging
                
def submit_16S_FAST(ar1,ar2,lr1,lr2,remark,name,result_path='/mnt/data/yangkai/16S_FAST/',group=0):
    result_path = os.path.join(result_path,remark,name)
    cmd_list = ['/root/anaconda3/bin/python3.7', 
                '/Bioinfo/16S-FAST-Tools/16S_FAST_Analysis_0.4.1.py', 
                '-ar1 ${ar1} -ar2 ${ar2} -lr1 ${lr1} -lr2 ${lr2}',
                '-o', '${outdir}',
                '-n', name,
                '-c', str(16), 
                '-minl', '1200', 
                '-maxl', '1700'
               ]
    if group == 1:
        cmd_list.append('-group')
    cmd=' '.join(cmd_list)
    config={ 
                'application':
                {
                    'cmd':cmd,
                    'id':1,
                    'name':'16S_FAST_Analysis_0.3.1',
                },
                'baseData':
                {
                    'remark':remark
                },
                'datas':
                [
                    {
                        'objectKey':ar1,'paramKey':'ar1'
                    },
                    {
                        'objectKey':ar2,'paramKey':'ar2'
                    },
                    {
                        'objectKey':lr1,'paramKey':'lr1'
                    },
                    {
                        'objectKey':lr2,'paramKey':'lr2'
                    },
                ],
                'name':name,
                'result':
                {
                    'path':result_path,
                    'paramKey':'outdir'
                },
                'type':'aln'    
    }
    config_json=(json.dumps(config))
    commond="curl -X POST https://bc.dev.germountx.com/api/tasks -H 'Content-Type: application/json' -d " + \
            "'" + config_json + "'"
    #logging.info(commond)
    out=os.popen(commond)
    logging.info('{} {}'.format(name,out.read()))
        
                
if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO) 
    date = time.strftime("%Y%m%d%H%M%S", time.localtime()) 
    __version__ = "0.1_20190918"
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version',version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('-ar1','--Assemble_Read_1',help='Assemble Read 1 fastq data absolute path in bc-input',required=True)
    parser.add_argument('-ar2','--Assemble_Read_2',help='Assemble Read 2 fastq data absolute path in bc-input',required=True)
    parser.add_argument('-lr1','--Linker_Read_1',help='Linker Read 1 fastq data absolute path in bc-input',required=True)
    parser.add_argument('-lr2','--Linker_Read_2',help='Linker Read 2 fastq data absolute path in bc-input',required=True)
    
    parser.add_argument('-n','--name',help='sample name.',required=True)
    parser.add_argument('-remark','--remark',help='batch name. default:Current Time',default=date)
    parser.add_argument('-o','--outdir',help='result dir in bc-output. default:/mnt/data/yangkai/16S_FAST/remark/name',\
                        default='/mnt/data/yangkai/16S_FAST/')
    args = parser.parse_args() 
    ar1 = args.Assemble_Read_1
    ar2 = args.Assemble_Read_2
    lr1 = args.Linker_Read_1
    lr2 = args.Linker_Read_2                
    remark = args.remark
    name = args.name
    
    submit_16S_FAST(ar1,ar2,lr1,lr2,remark,name)
    
