import os
import re
import argparse

def get_data_path_of_16S_FAST(data_info_file,path_out_File,Name_col,lib_col,P_col,L_col,data_path,pattern):
    f=data_info_file
    f_o = open(path_out_File,'w')
    total_info=open(f).readlines()
    info_col=[P_col,L_col]
    command ='qsctl ls qs://bc-input-prod' + data_path
    if pattern:
        command = command + ' | grep ' + pattern
    database=os.popen(command).read().split('\n')
    for sample in total_info:
        tmp = sample.split('\t')
        print(tmp[Name_col-1],end='\t',file=f_o)
        for col in info_col:
            match = ''
            if tmp[col-1]:
#                if tmp[lib_col-1]:
#                    match = '/' + tmp[lib_col-1] + '-'+ tmp[col-1]
#                else:
#                    match = '/' + tmp[col-1]
                match = '/' + tmp[col-1]
                for d in database:
                    if re.search(match,d):
                        print(d.split()[5] + '\t' + d.split()[3] + d.split()[4],end='\t',file=f_o)
            else:
                print('\t'*3,end='\t',file=f_o)
        print(file=f_o)              
    f_o.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-data_info_file','--data_info_file',help='data info file, Germount_16S_FAST_Data_Info.xlsx columns 7-11',required=True)
    parser.add_argument('-path_out_File','--path_out_File',help='out file',required=True)
    parser.add_argument('-Name_col','--Name_col',help='Sample Name columns number. default=2',type=int,default=2)
    parser.add_argument('-lib_col','--lib_col',help='DNA libary Name columns number. default=6',type=int,default=6)
    parser.add_argument('-P_col','--P_col',help='DNA Assemble libary Name columns number. default=7',type=int,default=7)
    parser.add_argument('-L_col','--L_col',help='DNA Linker libary Name columns number. default=11',type=int,default=11)
    parser.add_argument('-pattern','--pattern',help='mark of the batch in data(date,library-name...) . default=''',default='')
    parser.add_argument('-data_path','--data_path',help='data path in stoarge. default:/mnt/data/16S_FAST_data/',\
                        default='/mnt/data/16S_FAST_data/')
    args = parser.parse_args()
    get_data_path_of_16S_FAST(args.data_info_file,args.path_out_File,\
                              args.Name_col,args.lib_col,args.P_col,args.L_col,args.data_path,\
                             args.pattern)
 
