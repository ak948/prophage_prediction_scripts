import shutil

import pandas as pd
import os
import subprocess
import csv

df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid'].loc[10000:10350]


for id in biosampleid:
    subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)
    os.remove(f'{id}.tar.gz')
    os.chdir(f"{id}/phageboost_results")
    # Run prokka on predicted prophage genomes
    file_type = '.fna'
    base_name = f'{id}_p'
    file_num = 1
    for file in os.listdir():
        if os.path.isdir(f'{base_name}{file_num}'):
            os.chdir(f'{base_name}{file_num}')
            file_type = '.fna'
            base_name = f'{id}_p'
            if os.path.isfile(f'{base_name}{file_num}.fna') and not os.path.isfile(f'{base_name}{file_num}.gbk'):
                subprocess.run(
                    f"prokka --outdir {base_name}{file_num} --prefix {id}_p{file_num} --locus {id}_p{file_num} --cpus 40 --hmms ~/miniconda3/envs/akinsgenomad/db/hmm/all_phrogs.hmm {base_name}{file_num}{file_type} --compliant",
                    shell=True)
                os.chdir(f"{base_name}{file_num}")
                subprocess.run(
                    f"mv {id}_p{file_num}.gff {id}_p{file_num}.gbk {id}_p{file_num}.faa ..",
                    shell=True)
                os.chdir('..')
                subprocess.run(f" rm -r {base_name}{file_num}", shell=True)
                os.chdir('..')
                file_num += 1
            else:
                os.chdir('..')
                file_num += 1

    os.chdir('..')
    os.chdir('..')
    subprocess.run(f"tar -zcf {id}.tar.gz {id}", shell=True)
    shutil.rmtree(id)
