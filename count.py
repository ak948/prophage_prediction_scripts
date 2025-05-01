import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import pylab

df = pd.read_csv('prophage_results.tsv', sep='\t')
prophage_name = df['Prophage Name']

for name in prophage_name:
    with open('my_clusters.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t', fieldnames=['column1', 'column2'])
        for row in reader:
            if f'{name},' in row['column2']:
                cluster_rep = row['column1']
            elif f'{name}"' in row['column2']:
                cluster_rep = row['column1']
            elif row['column2'] == name:
                cluster_rep = row['column1']
    with open('Summary_taxonomy.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{cluster_rep}' in row['Genome']:
                cluster = row['Cluster']

    with open('Summary_taxonomy.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{name}' in row['Genome']:
                cluster_mapping[row['Genome']] = row['Cluster']
    if name in cluster_mapping:
        dt.at[df.index[dt['Genome'] == name], 'EPC Cluster'] = cluster_mapping[name]
