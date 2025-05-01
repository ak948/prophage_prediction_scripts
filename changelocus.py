import shutil

import pandas as pd
import os
import subprocess
import csv

os.chdir('Completed')
df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid'].loc[0:600]


for term in biosampleid:
    file_num = 1
    if os.path.isfile(f'{term}.tar.gz'):
        subprocess.run(f"tar -zxvf {term}.tar.gz {term}", shell=True)
        os.remove(f'{term}.tar.gz')
        os.chdir(f'{term}/phageboost_results')

        while True:
            folder_name = f"{term}_p{file_num}"
            folder_path = os.path.join(os.getcwd(), folder_name)

            if os.path.exists(folder_path):
                os.chdir(folder_path)
                file_name = f"{term}_p{file_num}"
                if os.path.isfile(f'{file_name}.fna'):
                    subprocess.run(
                        f"sed -i 's/gnl|Prokka|//gp' {term}_p{file_num}.fna",
                        shell=True)
                if os.path.isfile(f'{file_name}.gbk'):
                    subprocess.run(
                        f"sed -i 's/gnl|Prokka|//gp' {term}_p{file_num}.gbk",
                        shell=True)
                if os.path.isfile(f'{file_name}.gff'):
                    subprocess.run(
                        f"sed -i 's/>gnl|Prokka|//gp' {term}_p{file_num}.gff",
                        shell=True)

            if not os.path.exists(folder_path):
                break

            os.chdir('..')
            file_num += 1

        os.chdir('../..')
        subprocess.run(f"tar -zcvf {term}.tar.gz {term}", shell=True)
        shutil.rmtree(term)

