import pandas as pd
import argparse

def get_tax_from_mothur(tax_file,asv_file,tax_ratio_file,tax_reads_file):
    # get taxonmy from mothur result
    #in file
    #tax_file = '../2019_EQA/asv/merge_nonch.silva_132_v3.wang.taxonomy'
    #asv_file = '../2019_EQA/asv/asv_tab.txt'
    #out file
    #tax_ratio_file = '../2019_EQA/tax_ratio.xlsx'
    #tax_reads_file = '../2019_EQA/tax_reads.xlsx'
    
    df_count = pd.DataFrame()
    df_asv_reads = pd.read_csv(asv_file,sep='\t',index_col=0)
    df_asv_ratio = df_asv_reads / df_asv_reads.sum()
    df_count['Contigs']=df_asv_reads.sum()
    df_count['Asv']=(df_asv_ratio>0).sum()
    
    dict_tag = {'Kingdom':0,'Phylum':1,'Class':2,'Order':3,'Family':4,'Genus':5,'Species':6}
    df_asv_taxon = pd.read_csv(tax_file,sep='\t',header=None)
    df_asv_taxon.columns=(['Rep_Seq','Taxon'])
    for tax,pos in dict_tag.items():
        df_asv_taxon[tax] = df_asv_taxon['Taxon'].str.split(';').str[pos].str.split('(').str[0]
    df_asv_taxon.drop('Taxon',axis=1,inplace=True)
    df_asv_taxon['ID'] = pd.to_numeric(df_asv_taxon['Rep_Seq'].str.split('_').str[1])
    df_asv_taxon = df_asv_taxon.set_index('ID').drop('Rep_Seq',axis=1)
    
    df_taxon_ratio = pd.merge(df_asv_ratio,df_asv_taxon,left_index=True,right_index=True)
    df_taxon_reads = pd.merge(df_asv_reads,df_asv_taxon,left_index=True,right_index=True)
    
    writer_ratio = pd.ExcelWriter(tax_ratio_file)
    writer_reads = pd.ExcelWriter(tax_reads_file)
    
    
    for tax in dict_tag.keys():
        sample_name = list(df_asv_ratio.columns)
        sample_name.append(tax)
        df_tmp_ratio = df_taxon_ratio[sample_name].groupby([tax]).sum()
        df_tmp_reads = df_taxon_reads[sample_name].groupby([tax]).sum()
        unclass = df_tmp_ratio.index.get_level_values(level=0).str.contains(r'unclassified|other|unknown|uncultured|unidentified')
        df_tmp_ratio.loc['Unclassified']=df_tmp_ratio[unclass].sum()
        df_tmp_reads.loc['Unclassified']=df_tmp_reads[unclass].sum()
        unclass = df_tmp_ratio.index.get_level_values(level=0).str.contains(r'unclassified|other|unknown|uncultured|unidentified')
        df_tmp_ratio.drop(df_tmp_ratio[unclass].index,inplace=True)
        df_tmp_reads.drop(df_tmp_reads[unclass].index,inplace=True)
        
        df_tmp_ratio.sort_values(by=list(df_tmp_ratio.columns),ascending=False).to_excel(writer_ratio,tax)
        df_tmp_reads.sort_values(by=list(df_tmp_reads.columns),ascending=False).to_excel(writer_reads,tax)
        
        df_count[tax]=(df_tmp_ratio>0).sum()
    
    df_count.to_excel(writer_ratio,'Sta')
    df_count.to_excel(writer_reads,'Sta')
    writer_ratio.save()
    writer_reads.save()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-taxfile','--taxfile',help='mothur wang.taxonomy file',required=True)
    parser.add_argument('-table','--table',help='asv/otu table',required=True)
    parser.add_argument('-taxfile_reads','--taxfile_reads',\
                        help='output excel file for reads. default=tax_reads.xlsx',default='tax_reads.xlsx')
    parser.add_argument('-taxfile_ratio','--taxfile_ratio',\
                        help='output excel file for ratio. default=tax_ratio.xlsx',default='tax_ratio.xlsx')

    args = parser.parse_args()
    get_tax_from_mothur(args.taxfile,args.table,args.taxfile_ratio,args.taxfile_reads)
    
