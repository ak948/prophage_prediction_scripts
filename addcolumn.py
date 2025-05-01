import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import pylab
import gzip
from Bio import Entrez
import urllib
import time
import sys

csv.field_size_limit(sys.maxsize)

df = pd.read_csv('prophage_dataframe.tsv', sep='\t', dtype={'Maximum number of coding sequences in prophages in host': 'float', 'Average Prophage GC% in Host': 'float', 'Minimum Prophage GC% in Host': 'float'})
prophage_name = df['Prophage Name']

# cluster_reps = []
# cluster_qualitys = []
#
# for name in prophage_name:
#     # Get cluster for the prophage
#     # Get representative for the cluster and add to dictionary
#     # Get quality of the cluster representative and add to dictionary
#     # Make columns with the dictionarys as values
#     with open('prophage_dataframe.tsv', 'r') as file:
#         reader = csv.DictReader(file, delimiter='\t')
#         for row in reader:
#             if row['Prophage Name'] == name:
#                 cluster = row['Cluster']
#     with open('Summary_taxonomy.tsv', 'r') as file:
#         reader = csv.DictReader(file, delimiter='\t')
#         for row in reader:
#             if f'{cluster}' in row['Cluster']:
#                 cluster_rep = row['Genome']
#                 cluster_reps.append(cluster_rep)
#     with open('dbcheckv/quality_summary.tsv', 'r') as file:
#         reader = csv.DictReader(file, delimiter='\t')
#         for row in reader:
#             if f'{cluster_rep}' in row['contig_id']:
#                 cluster_quality = row['checkv_quality']
#                 cluster_qualitys.append(cluster_quality)
#
# df.insert(68, "Cluster Representative", cluster_reps)
# df.insert(69, "Cluster Quality", cluster_qualitys)

df.insert(68, "Cluster Representative", '')
df.to_csv('prophage_dataframe.tsv', sep='\t', index=False)
df.insert(69, "Cluster Quality", '')
df.to_csv('prophage_dataframe.tsv', sep='\t', index=False)

for name in prophage_name:
    with open('prophage_dataframe.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['Prophage Name'] == name:
                cluster = row['Cluster']
    with open('Summary_taxonomy.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{cluster}' in row['Cluster']:
                cluster_rep = row['Genome']
    with open('dbcheckv/quality_summary.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{cluster_rep}' in row['contig_id']:
                cluster_quality = row['checkv_quality']
    df.at[f'{name}', 'Cluster Representative'] = cluster_rep
    df.at[f'{name}', 'Cluster Quality'] = cluster_quality
    df.to_csv('prophage_dataframe.tsv', sep='\t', index=False)
