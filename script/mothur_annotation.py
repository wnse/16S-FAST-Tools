import os
import logging
db_path = '/Bioinfo/Database/'
cmd_path = '/Bioinfo/bin/'
mothur_db_fa = db_path + 'mothur/silva_132_v3.fa'
mothur_db_tax = db_path + 'mothur/silva_132_v3.2.tax'
mothur = cmd_path + 'mothur'
def mothur_annotation(infa,mothur_path,mothur_db,mothur_tax,threads=2):
    set_dir = os.dir
    cmd = mothur 
    cmd = cmd + ' "#classify.seqs(fasta=' + infa + ','
    cmd = cmd + 'template=' + mothur_db + ','
    cmd = cmd + 'taxonomy=' + mothur_tax + ','
    cmd = cmd + 'processors=' + str(threads) + ')"'
    os.chdir(os.path.dirname(infa))
    status=os.system(cmd)
    if status==0:
        logging.info(' done')
    else:
        logging.info(' exit')
        
if __name__ == '__main__':
    infa = '../test/merged.asv.fasta'
    threads = 4
    status =  mothur_annotation(infa,mothur,mothur_db_fa,mothur_db_tax,threads)
    