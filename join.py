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

# Make new file to merge into, change headers as needed
with open('datatoadd.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['Prophage Name', 'Cluster Representative', 'Cluster Quality', 'Vcontact Cluster',
                         'Genomes in VC', 'Representative Genomes in VC', 'Total Genomes in VC'])

for name in prophage_name:
    with open('prophage_dataframe.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['Prophage Name'] == name:
                cluster = row['Cluster']
                cluster_rep = row['Cluster Representative']
    # with open('Summary_taxonomy.tsv', 'r') as file:
    #     reader = csv.DictReader(file, delimiter='\t')
    #     for row in reader:
    #         if f'{cluster}' in row['Cluster']:
    #             cluster_rep = row['Genome']
    with open('dbcheckv/quality_summary.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{cluster_rep}' in row['contig_id']:
                cluster_quality = row['checkv_quality']
    with open('vcontact_clusters.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{cluster},' in row['EPCs Present']:
                vc = row['VC']
                vc_genomes = row['Genomes in VC']
                vc_representatives = row['Representative Genomes in VC']
                vc_prophages = row['Total Prophage Genomes in VC']
                vc_totals = int(vc_prophages) + int(vc_genomes) - int(vc_representatives)
    # Add collected data to new tsv file
    with open('datatoadd.tsv', 'a') as outfile:
        tsv_writer = csv.writer(outfile, delimiter='\t')
        tsv_writer.writerow([f'{name}', f'{cluster_rep}', f'{cluster_quality}', f'{vc}', f'{vc_genomes}',
                             f'{vc_representatives}', f'{vc_totals}'])

# #Open main file as df
# df = pd.read_csv('prophage_dataframe.tsv', sep='\t', low_memory=False)
# # Open new tsv file as dt
# dt = pd.read_csv('datatoadd.tsv', sep='\t')
# # Add columns in new file to main file
# dg = df.join(dt.set_index('Prophage Name'), on='Prophage Name')
# dg.to_csv('prophage_dataframe.tsv', sep='\t', index=False)

# import pandas as pd
# dt = pd.read_csv('removed.txt', sep='\t')
# biosample_ids = dt['BiosampleID']
# df = pd.read_csv('amr_genes.tsv', sep='\t')
# dg = df[df.name.isin(biosample_ids) == False]
# dg.to_csv('amr_genes.tsv', sep='\t', index=False)
