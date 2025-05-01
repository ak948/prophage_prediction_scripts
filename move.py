import shutil
import pandas as pd
import os
import subprocess
import csv

df = pd.read_excel('genomeinfo.tsv', sheet_name='prokaryotes', header=0)
biosampleid = df['BioSample'].loc[0:500]
for id in biosampleid:
    subprocess.run(f"cp {id}.tar.gz ~/Documents/copy2laptop", shell=True)