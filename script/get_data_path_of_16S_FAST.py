import os
import re
import argparse


def get_data_path_of_16S_FAST(data_info_file,
                              path_out_File,
                              Name_col,
                              P_col,
                              L_col,
                              data_path,
                              pattern):
    '''从对象存储获取数据路径，输入文件中包含样本名称与A库和L库的数据名称
    参数：
        data_info_file: 输入文件，包含样本名称、A库数据名称、L库数据名称
        path_out_File: 输出文件
        Name_col: 样本名称列
        P_col: A库数据名称
        L_col: L库数据名称
        data_path: 对象存储数据路径
        pattern: 匹配标识（以防在存储中包含相同的数据名称，指定的标识，可以使用数据上传日期作为标识）
    返回：
        无
    '''
    f = data_info_file
    f_o = open(path_out_File, 'w')
    total_info = open(f).readlines()
    info_col = [P_col, L_col]
    command = 'qsctl ls qs://bc-input-prod' + data_path
    if pattern:
        command = command + ' | grep ' + pattern
    database = os.popen(command).read().split('\n')
    for sample in total_info:
        tmp = sample.split('\t')
        print(tmp[Name_col - 1], end='\t', file=f_o)
        for col in info_col:
            if tmp[col - 1]:
                match = '/' + tmp[col - 1]
                for d in database:
                    if re.search(match, d):
                        print(d.split()[5] + '\t' + d.split()[3] + d.split()[4], end='\t', file=f_o)
            else:
                print('\t' * 3, end='\t', file=f_o)
        print(file=f_o)
    f_o.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--data_info_file',
                        help='data info file',
                        required=True)
    parser.add_argument('-o', '--path_out_File',
                        help='out file',
                        required=True)
    parser.add_argument('-n', '--Name_col',
                        help='Sample Name columns number.',
                        type=int,
                        default=2)
    # parser.add_argument('-lib_col',
    #                     '--lib_col', help='DNA libary Name columns number.',
    #                     type=int,
    #                     default=6)
    parser.add_argument('-P', '--P_col',
                        help='DNA Assemble libary Name columns number.',
                        type=int,
                        default=9)
    parser.add_argument('-L', '--L_col', help='DNA Linker libary Name columns number.',
                        type=int,
                        default=13)
    parser.add_argument('-p', '--pattern',
                        help='mark of the batch in data(date).',
                        default='')
    parser.add_argument('-data_path', '--data_path',
                        help='data path in stoarge.',
                        default='/mnt/data/16S_FAST_data/')
    args = parser.parse_args()
    get_data_path_of_16S_FAST(args.data_info_file,
                              args.path_out_File,
                              args.Name_col,
                              args.P_col,
                              args.L_col,
                              args.data_path,
                              args.pattern)
