import os
import csv
import pandas as pd

df = pd.read_csv('amr_genes.tsv', sep='\t')
unique_AMR = df['AMR Gene Family'].unique()
unique_AMR = sorted(unique_AMR)  # Sort unique values for consistent column order

dt = pd.read_csv('virulence_genes.tsv', sep='\t')
unique_VIR = dt['Gene'].unique()
unique_VIR = sorted(unique_VIR)  # Sort unique values for consistent column order

with open('argvirtypes.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    header = ["BiosampleID", "Prophage AMR Gene Count"]
    header.extend([f"Prophage ARG Family '{arg}' Count" for arg in unique_AMR])
    header.extend(["Prophage Virulence Gene Count"])
    header.extend([f"Prophage Virulence Gene '{virulence}' Count" for virulence in unique_VIR])
    tsv_writer.writerow(header)


df = pd.read_csv('results.tsv', sep='\t')
genome_name = df['BiosampleID']

for name in genome_name:
    amr_dict = {}
    vir_dict = {}
    for amr in unique_AMR:
        amr_dict[amr] = 0
    for vir in unique_VIR:
        vir_dict[vir] = 0
    amr_count = 0
    vir_count = 0
    with open('amr_genes.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{name}_p' in row['Prophage Name']:
                amr_count += 1
                amr_type = row['AMR Gene Family']
                amr_dict[amr_type] += 1

    with open('virulence_genes.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{name}_p' in row['Prophage ID']:
                vir_count += 1
                vir_type = row['Gene']
                vir_dict[vir_type] += 1


    with open('argvirtypes.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        row = [f'{name}', f'{amr_count}']
        row.extend([f'{amr_dict[arg]}' for arg in unique_AMR])
        row.extend([f'{vir_count}'])
        row.extend([f'{vir_dict[vir]}' for vir in unique_VIR])
        tsv_writer.writerow(row)
