import os
import pandas as pd
import re
import argparse
def merge_tax(indir, filename='tax.ratio.xlsx'):
    s_dir = [os.path.join(indir, i) for i in os.listdir(indir)]
    df_dict = {}
    for s in s_dir:
        t_dir = os.path.join(s, 'result')
        for file in os.listdir(t_dir):
            if re.search(filename, file):
                f = os.path.join(t_dir, file)
                for tax in dict_tag.keys():
                    df_dict[tax] = df_dict.get(tax, pd.DataFrame())
                    df_tmp = pd.read_excel(f, tax, index_col=0).drop(tax, axis=1)
                    if df_dict[tax].empty:
                        df_dict[tax] = df_tmp
                    else:
                        df_dict[tax] = pd.merge(df_dict[tax], 
                                                df_tmp, 
                                                left_on='Taxonomy', 
                                                right_on='Taxonomy', 
                                                how='outer')
                df_dict['Sta'] = df_dict.get('Sta', pd.DataFrame())
                df_tmp = pd.read_excel(f, 'Sta', index_col=0).rename(index={'test':s})
                df_dict['Sta'] = pd.concat([df_dict['Sta'], df_tmp])
    return df_dict

def write_file(df_dict, outfile):
    writer_ratio = pd.ExcelWriter(outfile)
    for tax in dict_tag.keys():
        df_tmp = df_dict[tax].fillna(0)
        df_tmp.insert(0,tax,df_tmp['Taxonomy'].str.split('__').str[-1])
        df_tmp.set_index(tax).sort_index().to_excel(writer_ratio, tax)
    df_dict['Sta'].to_excel(writer_ratio, 'Sta')
    writer_ratio.save()
                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--indir', 
                        help='result dir which including all sample result dir',
                        required=True)
    parser.add_argument('-or', '--outRatiofile', 
                        help='output excel file of total statistics info',
                        default='total_tax_ratio.xlsx')
    parser.add_argument('-oc', '--outCountfile', 
                        help='output excel file of total statistics info',
                        default='total_tax_count.xlsx')
    args = parser.parse_args()
    dict_tag = {'Kingdom':0,'Phylum':1,'Class':2,'Order':3,'Family':4,'Genus':5,'Species':6}
    df_dict = merge_tax(args.indir, filename='tax.ratio.xlsx')
    write_file(df_dict, args.outRatiofile)
    df_dict = merge_tax(args.indir, filename='tax.reads.xlsx')
    write_file(df_dict, args.outCountfile)
