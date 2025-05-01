import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np

with open('bacteriaepcs.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(["BiosampleID", "Phage EPCs in Genome"])

df = pd.read_csv('results.tsv', sep='\t')
biosampleid = df['BiosampleID']

for id in biosampleid:
    epc_dict = []

    with open('prophage_dataframe.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['BiosampleID'] == id:
                epc_dict.append(row['Cluster'])

    with open('bacteriaepcs.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow([f'{id}', f'{epc_dict}'])

    print(f'{id} Done')