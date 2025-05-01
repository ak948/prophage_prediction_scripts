import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import pylab
import re


with open('vcontact_clusters.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['VC', 'Genomes in VC', 'Representative Genomes in VC', 'Total Prophage Genomes in VC',
                         'Number of Overlapped Genomes', 'EPCs Present'])

df = pd.read_csv('vcontact_genomes.csv', sep=',',)
genome_name = df['Genome']

row_number = 0

for name in genome_name:

    vc_value = []
    epc_genomes = 0
    rep_genome = 0
    epc_clusters = []
    vcs = []
    subcluster_value = []
    overlapped_genomes = 0
    vc_genome = 0
    row_number += 1

    with open('vcontact_genomes.csv', 'r') as file:
        reader = csv.DictReader(file, delimiter=',')

        for row in reader:
            if name == row['Genome']:
                if row['VC Status'] == "Clustered":
                    # Get VC value
                    vc_value = 'VC_' + row['VC']
                    subcluster_value = row['VC Subcluster']
                    vc_genome += 1
                    if "SAM" in row['Genome']:
                        epc_genomes = row['Genomes in EPC']
                        rep_genome += 1
                        epc_clusters = row['EPC Cluster']
                elif row['VC Status'] == "Outlier":
                    vc_status = row['VC Status']
                    vc_value = f'{name}_{vc_status}'
                    vc_genome += 1
                    if "SAM" in row['Genome']:
                        epc_genomes = row['Genomes in EPC']
                        rep_genome += 1
                        epc_clusters = row['EPC Cluster']
                elif row['VC Status'] == "Singleton":
                    vc_status = row['VC Status']
                    vc_value = f'{name}_{vc_status}'
                    vc_genome += 1
                    if "SAM" in row['Genome']:
                        epc_genomes = row['Genomes in EPC']
                        rep_genome += 1
                        epc_clusters = row['EPC Cluster']
                elif row['VC Status'] == "Clustered/Singleton":
                    # Get VC value
                    vc_value = 'VC_' + row['VC']
                    subcluster_value = row['VC Subcluster']
                    vc_genome += 1
                    if "SAM" in row['Genome']:
                        epc_genomes = row['Genomes in EPC']
                        rep_genome += 1
                        epc_clusters = row['EPC Cluster']
                elif "Overlap" in row['VC Status']:
                    vc_value = f'{name}_overlap'
                    vc_genome += 1
                    if "SAM" in row['Genome']:
                        epc_genomes = row['Genomes in EPC']
                        rep_genome += 1
                        epc_clusters = row['EPC Cluster']

        dt = pd.read_csv('vcontact_clusters.tsv', sep='\t')
        if vc_value in dt['VC'].values:
            for index, row in dt.iterrows():
                if row['VC'] == vc_value:
                    index_row = index
            vc_genomes = dt.at[index_row, 'Genomes in VC']
            vc_reps = dt.at[index_row, 'Representative Genomes in VC']
            vc_props = dt.at[index_row, 'Total Prophage Genomes in VC']
            vc_epcs = dt.at[index_row, 'EPCs Present']
            dt.at[index_row, 'Genomes in VC'] = 1 + int(vc_genomes)
            dt.at[index_row, 'Representative Genomes in VC'] = int(rep_genome) + int(vc_reps)
            dt.at[index_row, 'Total Prophage Genomes in VC'] = int(epc_genomes) + int(vc_props)
            if len(epc_clusters) > 2:
                dt.at[index_row, 'EPCs Present'] = vc_epcs + epc_clusters + ','
            dt.to_csv('vcontact_clusters.tsv', sep='\t', index=False)
        else:
            with open('vcontact_clusters.tsv', 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow([f'{vc_value}', f'{vc_genome}', f'{rep_genome}', f'{epc_genomes}', 0, f'{epc_clusters},'])

