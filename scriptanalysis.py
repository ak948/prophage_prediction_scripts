import shutil

import pandas as pd
import os
import subprocess
import csv


# Get list of biosampleIDs
df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid'].loc[0:10]
# Write headings for file results for each genome will be recorded in
with open('results.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(
        ['Biosampleid', 'Species', 'Collection Date', 'Collection Location', 'Isolation Source', 'Number of Prophages',
         'Average Prophage Length', 'Minimum Prophage Length', 'Maximum Prophage Length',
         'Average number of prophage coding sequences', 'Minimum number of prophage coding sequences',
         'Maximum number of prophage coding sequences', 'Number of prophages containing AMR Genes',
         'Total number of AMR Genes in prophages'])
# Write headings for file the rgi-card results will be recorded in for each genome
with open('amr_genes.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['AMR Gene ID', 'ORF ID', 'Source Sequence', 'Cut-Off', 'Drug Class', 'Resistance Mechanism',
                         'AMR Gene Family', '% Length of Reference Genome'])

total_lengths = []
total_genes = []
total_amr_genes = []

from analysis import *

for id in biosampleid:
    subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)
    extractinfo = extract_info(id) # Get info on the genome from the results.tsv table
    predictioninfo = prediction_info(id) # Get info on the predicted prophages
    amrgenes = amr_genes(id) # Use RGI-CARD to predict AMR genes
    #closestrelatives
    #prophageclustering
    #comparativegenomics
    addinfo = add_info(id) # Add info on genome and prophages to table
    subprocess.run(f"tar -zcf {id}.tar.gz {id}", shell=True)
    shutil.rmtree(id)

creategraphs

