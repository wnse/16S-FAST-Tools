import pandas as pd
import numpy as np
import os
import sys
import time
import logging
import argparse
import submit_16S_FAST


def submit_16S_FAST_group_from_list(df_tmp, colums_tmp, remark_tmp, name_tmp, group=1):
    for i in df_tmp.index:
        sample_name = df_tmp.loc[i, int(name_tmp) - 1]
        submit_16S_FAST.submit_16S_FAST(df_tmp.loc[i, int(colums_tmp[0]) - 1],
                                        df_tmp.loc[i, int(colums_tmp[1]) - 1],
                                        df_tmp.loc[i, int(colums_tmp[2]) - 1],
                                        df_tmp.loc[i, int(colums_tmp[3]) - 1],
                                        remark_tmp,
                                        sample_name,
                                        group=group)
        logging.info('{}'.format(sample_name))


if __name__ == '__main__':
    date = time.strftime("%Y%m%d%H%M%S", time.localtime())
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-file', '--list_file', help='16S_data_list', required=True)
    parser.add_argument('-data_col', '--col_of_data_path', help='columns of data path (a1,a2,l1,l2).',
                        default='2,4,6,8')
    parser.add_argument('-name_col', '--col_of_name', help='columns of sample name.', default=1)
    parser.add_argument('-remark', '--remark', help='batch name. default:Current Time', default=date)
    parser.add_argument('-group', '--group', help='keep group umi seq in result', default=1, type=int, choices=([0, 1]))
    args = parser.parse_args()

    df = pd.read_csv(args.list_file, header=None, sep='\t')
    colums = args.col_of_data_path.split(',')
    remark = args.remark
    name = args.col_of_name
    group = args.group

    submit_16S_FAST_group_from_list(df, colums, remark, name, group=group)
