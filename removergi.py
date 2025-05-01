import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil

# os.chdir('Completed')

#Get list of biosampleIDs
df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid'].loc[0:100]

for id in biosampleid:

    if os.path.isfile(f'{id}.tar.gz'):
        subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)
        os.remove(f'{id}.tar.gz')
        os.chdir(f'{id}/phageboost_results')
        prophage_num = 0

        for filename in os.listdir():
            if filename.startswith(f'{id}_p'):
                prophage_num += 1
                os.chdir(f'{id}_p{prophage_num}')

                for filename in os.listdir():
                    if filename.endswith('rgi.json'):
                        os.remove(f'{id}_p{prophage_num}_rgi.json')
                    if filename.endswith('rgi.txt'):
                        os.remove(f'{id}_p{prophage_num}_rgi.txt')
                    keyword = 'tmp'
                    if keyword in filename:
                        subprocess.run('rm *tmp*', shell=True)
                    if filename.endswith('.temp'):
                        subprocess.run('rm *.temp* *tmp*', shell=True)
                    if filename.endswith('output.tsv'):
                        os.remove('output.tsv')
                    if filename.endswith('test.db'):
                        os.remove('test.db')
                os.chdir('..')



    os.chdir('/home/ggb_andrew/Documents/prophage_prediction')
    subprocess.run(f"tar -zcvf {id}.tar.gz {id}", shell=True)
    shutil.rmtree(id)