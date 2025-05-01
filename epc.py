import csv
import pandas as pd
import sys
csv.field_size_limit(sys.maxsize)

with open('datatoadd.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Prophage Name', 'Cluster', 'Cluster Representative'])

df = pd.read_csv('prophage_dataframe.tsv', sep='\t')
prophage_name = df['Prophage Name']


for name in prophage_name:
    cluster_rep = None
    votu = None
    with open('new_my_clusters.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{name},' in row['Genomes']:
                cluster_rep = row['Representative']
            elif f'{name}n' in row['Genomes']:
                cluster_rep = row['Representative']
            elif row['Genomes'] == name:
                cluster_rep = row['Representative']
    with open('votu.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['Representative'] == cluster_rep:
                votu = row['vOTU']

    with open('datatoadd.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow([f'{name}', f'{votu}', f'{cluster_rep}'])

        # Cluster Quality, Species/Genus, Vcontact stuff,


