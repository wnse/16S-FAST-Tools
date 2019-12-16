import os
import pandas as pd
import re
import numpy as np
import argparse
def merge_sta(indir):
    s_dir = [os.path.join(indir, i) for i in os.listdir(indir)]
    df = pd.DataFrame()
    for s in s_dir:
        t_dir = os.path.join(s, 'result')
        for file in os.listdir(t_dir):
            if re.search('sta\.txt', file):
                f = os.path.join(t_dir, file)
                if df.empty:
                    df = pd.read_csv(f, sep='\t')
                else:
                    df = pd.merge(df, 
                                  pd.read_csv(f, 
                                              sep='\t', 
                                              index_col=0).rename(columns={'test':s}),
                                  left_on='description',
                                  right_on='description',
                                  how='outer')
    df['ID'] = np.array(df['description'].str.split().str[0]).astype(np.int32)
    return df.set_index('ID').sort_index()
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--indir', 
                        help='result dir which including all sample result dir',
                        required=True)
    parser.add_argument('-o', '--outfile', 
                        help='output excel file of total statistics info',
                        default='total_sta.xlsx')
    args = parser.parse_args()
    df = merge_sta(args.indir)
    df.to_excel(args.outfile, index=False)
