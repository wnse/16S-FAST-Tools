import argparse
import pandas as pd
from Bio import SeqIO

def get_consensus_seq_from_mothur(taxfile, fafile):
    df = pd.read_csv(taxfile, sep='\t', header=None)
    df.columns=['ID','Tax']
    df['umiID'] = df['ID'].str.extract(r'(\S+)_\d+')
    df['Tax'] = df['Tax'].replace('\(\d+\)','',regex=True)
    
    seq_len = {}
    for record in SeqIO.parse(fafile,'fasta'):
        seq_len[str(record.id)] = len(record.seq)
    df_fa_len = pd.DataFrame.from_dict(seq_len,orient='index',
                                       columns=(['Seq_Len']))
    df = pd.merge(df,df_fa_len,left_on='ID',right_index=True)
    
    consensus_umi = {}
    unconsensus_umi = {}
    for i in df[['umiID','Tax']].to_dict('split')['data']:
        consensus_umi[i[0]] = consensus_umi.get(i[0],[])
        consensus_umi[i[0]].append(i[1])
    for i,v in consensus_umi.items():
        if len(set(v)) > 1:
            unconsensus_umi[i] = v
    for i,v in unconsensus_umi.items():
        del consensus_umi[i]
    df_unconsensus = df.set_index('umiID').loc[unconsensus_umi.keys()]
    df_unconsensus = df_unconsensus[['ID','Seq_Len','Tax']]
    df_unconsensus = df_unconsensus.sort_values(by='ID')
    df_consensus = df[df['umiID'].isin(consensus_umi.keys())].groupby(
        ['umiID']).apply(lambda t: t[t.Seq_Len==t.Seq_Len.max()])
    df_consensus = df_consensus[['ID','Seq_Len','Tax']].set_index('ID')

    return df_consensus, df_unconsensus

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tax','--taxonomy',required=True,
                       help='mothur taxonomy file')
    parser.add_argument('-fa','--fasta',required=True,
                       help='fasta file used to taxonomy')
    args = parser.parse_args()
    df1, df2 = get_consensus_seq_from_mothur(args.taxonomy, args.fasta)
    df1.to_csv(args.taxonomy+'.consensus.csv',sep='\t')
    df2.to_csv(args.taxonomy+'.unconsensus.csv',sep='\t')
