# Group ARGs by Drug Class, filter by ;
# Group Phage-Defence systems by Defense Type
# Group virulence gene by name
# Paired T-test for overall ARGs, and for each drug class, as bacteria count is dependent on prophage count?

# Make lists of Drug Classes from bacteria_amr_genes, put each AMR Gene Family in the lists for each of its Drug Class mentioned
# For each Drug Class have dictionary for prophage and bacteria
# In results.tsv go row by row, for each AMR Gene Family in Drug Class, look for "Bacteria ARG Family '{gene}' Count" and same for prophage, adding them to drug class count
# Put in output tsv with columns for prophage and bacteria count for each drug class, and rows for each biosampleid
# Go onto next Drug Class, once done with all Drug Classes, Defence Types, Virulence Genes, go onto next row

import os
import pandas as pd
import csv

dt = pd.read_csv('genes_by_drug_class.tsv', sep='\t')
dg = pd.read_csv('defenses_by_defense_type.tsv', sep='\t')
dh = pd.read_csv('../../Count/bacteria_virulence_genes.tsv', sep='\t')
virulence_list = dh['Gene'].unique()

with open('prophage_gene_impact_output.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    header = ["BiosampleID", "Bacteria ARG Count", "Prophage ARG Count", "Bacteria Defence System Count", "Prophage Defence System Count", "Bacteria Virulence Gene Count", "Prophage Virulence Gene Count"]
    header.extend([f"Bacteria ARG '{column}' Count" for column in dt.columns])
    header.extend([f"Prophage ARG '{column}' Count" for column in dt.columns])
    header.extend([f"Bacteria Defence Type'{column}' Count" for column in dg.columns])
    header.extend([f"Prophage Defence Type '{column}' Count" for column in dg.columns])
    header.extend([f"Bacteria Virulence Gene '{column}' Count" for column in virulence_list])
    header.extend([f"Prophage Virulence Gene '{column}' Count" for column in virulence_list])
    tsv_writer.writerow(header)

bacteria_arg_count = 0
prophage_arg_count = 0
bacteria_defence_count = 0
prophage_defence_count = 0
bacteria_virulence_count = 0
prophage_virulence_count = 0
row_num = 0
with open('../results.tsv', 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        df = pd.read_csv('../results.tsv', sep='\t')
        biosampleid = df['BiosampleID']
        genome_id = biosampleid[row_num]
        prophage_arg_count = df[f"Prophage AMR Gene Count"][row_num]
        bacteria_arg_count = df[f"Bacteria ARG Count"][row_num]
        bacteria_defence_count = df[f"Bacteria Phage Defense System Count"][row_num]
        prophage_defence_count = df[f"Prophage Defense System Count"][row_num]
        bacteria_virulence_count = df[f"Bacteria Virulence Gene Count"][row_num]
        prophage_virulence_count = df[f"Prophage Virulence Gene Count"][row_num]
        dt = pd.read_csv('genes_by_drug_class.tsv', sep='\t')
        arg_columns = dt.columns
        bacteria_arg_counts = {}
        prophage_arg_counts = {}
        for column in arg_columns:
            bacteria_arg_counts[column] = 0
            prophage_arg_counts[column] = 0
            column_genes = dt[column]
            for gene in column_genes:
                if pd.isna(gene) or gene == '':
                    continue
                else:
                    bacteria_arg_counts[column] += df[f"Bacteria ARG Family '{gene}' Count"][row_num]
                    column_name = f"Prophage ARG Family '{gene}' Count"
                    if column_name in df.columns:
                        prophage_arg_counts[column] += df[f"Prophage ARG Family '{gene}' Count"][row_num]

        dg = pd.read_csv('defenses_by_defense_type.tsv', sep='\t')
        defence_columns = dg.columns
        bacteria_defence_counts = {}
        prophage_defence_counts = {}
        for column in defence_columns:
            bacteria_defence_counts[column] = 0
            prophage_defence_counts[column] = 0
            column_genes = dg[column]
            for gene in column_genes:
                if pd.isna(gene) or gene == '':
                    continue
                else:
                    bacteria_defence_counts[column] += df[f"Bacteria Phage Defense System Subtype '{gene}' Count"][row_num]
                    column_name = f"Prophage Phage Defense System Subtype '{gene}' Count"
                    if column_name in df.columns:
                        prophage_defence_counts[column] += df[f"Prophage Phage Defense System Subtype '{gene}' Count"][row_num]

        dh = pd.read_csv('../../Count/bacteria_virulence_genes.tsv', sep='\t')
        virulence_list = dh['Gene'].unique()
        bacteria_virulence_counts = {}
        prophage_virulence_counts = {}
        for column in virulence_list:
            bacteria_virulence_counts[column] = 0
            prophage_virulence_counts[column] = 0
            bacteria_virulence_counts[column] += df[f"Bacteria Virulence Gene '{column}' Count"][
                row_num]
            column_name = f"Prophage Virulence Gene '{column}' Count"
            if column_name in df.columns:
                prophage_virulence_counts[column] += df[f"Prophage Virulence Gene '{column}' Count"][
                    row_num]

        row_num +=1

        with open('prophage_gene_impact_output.tsv', 'a') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            output_row = [f'{genome_id}', f'{bacteria_arg_count}', f'{prophage_arg_count}', f'{bacteria_defence_count}', f'{prophage_defence_count}', f'{bacteria_virulence_count}', f'{prophage_virulence_count}']
            output_row.extend([f'{bacteria_arg_counts[column]}' for column in arg_columns])
            output_row.extend([f'{prophage_arg_counts[column]}' for column in arg_columns])
            output_row.extend([f'{bacteria_defence_counts[column]}' for column in defence_columns])
            output_row.extend([f'{prophage_defence_counts[column]}' for column in defence_columns])
            output_row.extend([f'{bacteria_virulence_counts[column]}' for column in virulence_list])
            output_row.extend([f'{prophage_virulence_counts[column]}' for column in virulence_list])
            tsv_writer.writerow(output_row)



