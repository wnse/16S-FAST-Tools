from Bio import SeqIO
import pandas as pd
import argparse

def get_asv_seq_from_fasta(rawf,name,outf,tag):
    asv_fasta = open(outf,'w')
    seq_size = {}
    seq_id = {}
    n=0
    for seq in SeqIO.parse(rawf,'fasta'):
        s = str(seq.seq)
        seq_id[s] = seq_id.get(s,[])
        seq_id[s].append(str(seq.id))
    tmp=[]
    for seq,seqid in sorted(seq_id.items(),key = lambda d:len(d[1]),reverse=True):
        size = len(seqid)
        print ('>' + name + '_' + str(n) + ' ' + tag + ';size=' + str(size),file=asv_fasta)
        print (seq,file=asv_fasta)
        tmp.extend([[s,n,size] for s in seq_id[seq]])
        n+=1
    asv_fasta.close()
    df = pd.DataFrame(tmp)
    df.columns=(['Seq_ID','ASV_ID','ASV_Size'])
    return df
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fa','--fasta',help='input fasta file',required=True)
    parser.add_argument('-name','--name',help='sample name',required=True)
    parser.add_argument('-outfa','--out_fasta',\
                        help='output fasta file. default=name.asv.fa',default='')
    parser.add_argument('-tag','--tag',\
                        help='tag in fasta file. default=asv',default='asv')

    args = parser.parse_args()
    out_file = ''
    if args.out_fasta:
        out_file = args.out_fasta
    else:
        out_file = args.name + '.asv.fa'
    df = get_asv_seq_from_fasta(args.fasta,args.name,out_file,args.tag)
    df.to_csv(out_file+'.idinfo',sep='\t',index=False)
    
