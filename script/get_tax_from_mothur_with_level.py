import pandas as pd
import argparse

def get_tax_from_mothur_with_level(tax_file,asv_file,tax_ratio_file,tax_reads_file):
    df_count = pd.DataFrame()
    df_asv_reads = pd.read_csv(asv_file,sep='\t',index_col=0)
    df_asv_ratio = df_asv_reads / df_asv_reads.sum()
    df_count['Contigs']=df_asv_reads.sum()
    df_count['Asv']=(df_asv_ratio>0).sum()
    dict_tag = {'Kingdom':0,'Phylum':1,'Class':2,'Order':3,'Family':4,'Genus':5,'Species':6}
    dict_tag_r = dict(list((v,k) for k,v in dict_tag.items()))
    df_asv_taxon = pd.read_csv(tax_file,sep='\t',header=None, index_col=0)
    df_asv_taxon.columns=(['Taxon'])
    for pos,tax in dict_tag_r.items():
        tmp_list=list(tax[0].lower() for i in range(df_asv_taxon.shape[0]))
        df_asv_taxon[tax] = tmp_list
        df_asv_taxon[tax] = df_asv_taxon[tax].str.cat(df_asv_taxon['Taxon'].str.replace('\(\d+\)','').str.split(';').str[pos],sep='__')
        if pos > 0:
            df_asv_taxon[tax] = df_asv_taxon[dict_tag_r[pos-1]].str.cat(df_asv_taxon[tax],sep=';')
    df_asv_taxon.drop('Taxon',axis=1,inplace=True)
    #df_asv_taxon['ID'] = pd.to_numeric(df_asv_taxon['Rep_Seq'].str.split('_').str[1])
    #df_asv_taxon = df_asv_taxon.set_index('ID').drop('Rep_Seq',axis=1)
    df_taxon_ratio = pd.merge(df_asv_ratio,df_asv_taxon,left_index=True,right_index=True)
    df_taxon_reads = pd.merge(df_asv_reads,df_asv_taxon,left_index=True,right_index=True)
    writer_ratio = pd.ExcelWriter(tax_ratio_file)
    writer_reads = pd.ExcelWriter(tax_reads_file)
    sample_name = list(df_asv_ratio.columns)
    for tax in dict_tag.keys():
        sample_name_tmp = list(df_asv_ratio.columns)
        sample_name_tmp.append(tax)
        df_tmp_ratio = df_taxon_ratio[sample_name_tmp].groupby([tax]).sum().reset_index().rename(columns={tax:'Taxonomy'})
        df_tmp_ratio.insert(0,tax,df_tmp_ratio['Taxonomy'].str.split('__').str[-1])
        df_tmp_reads = df_taxon_reads[sample_name_tmp].groupby([tax]).sum().reset_index().rename(columns={tax:'Taxonomy'})
        df_tmp_reads.insert(0,tax,df_tmp_ratio['Taxonomy'].str.split('__').str[-1])
        df_tmp_ratio.sort_values(by=sample_name,ascending=False).reset_index().drop('index',axis=1).to_excel(writer_ratio,tax)
        df_tmp_reads.sort_values(by=sample_name,ascending=False).reset_index().drop('index',axis=1).to_excel(writer_reads,tax)
        df_count[tax]=(df_tmp_ratio[sample_name]>0).sum()
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
    get_tax_from_mothur_with_level(args.taxfile,args.table,args.taxfile_ratio,args.taxfile_reads)
