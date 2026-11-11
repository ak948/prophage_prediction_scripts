# Look at each gene in prophages and see if there are copies of the same gene elsewhere in the bacteria too
# Do for each gene separately, and also by grouping ARGs by drug resistances, and defence systems by their type
# For each category, say how many isolates have no copies of the gene, how many have only in bacteria, of ones with copy in prophage,
# how many have 0 additional copies, 1 additional copy, 2 copies, 3+
# Do for each gene/drug class separately

import os
import pandas as pd
import csv

dt = pd.read_csv('genes_by_drug_class.tsv', sep='\t')
dg = pd.read_csv('defenses_by_defense_type.tsv', sep='\t')
dh = pd.read_csv('bacteria_virulence_genes.tsv', sep='\t')
virulence_list = dh['Gene'].unique()
di = pd.read_csv('bacteria_amr_genes.tsv', sep='\t')
arg_family_list = di['AMR Gene Family'].unique()
dj = pd.read_csv('bacteria_phage_defenses.tsv', sep='\t')
defence_subtype_list = dj['Defense Subtype'].unique()

with open('prophage_gene_impact_output.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    header = ["BiosampleID", "Bacteria ARG Count", "Prophage ARG Count", "Bacteria Defence System Count", "Prophage Defence System Count", "Bacteria Virulence Gene Count", "Prophage Virulence Gene Count"]
    header.extend([f"Bacteria ARG Drug Class '{column}' Count" for column in dt.columns])
    header.extend([f"Prophage ARG Drug Class '{column}' Count" for column in dt.columns])
    header.extend([f"Bacteria ARG '{column}' Count" for column in arg_family_list])
    header.extend([f"Prophage ARG '{column}' Count" for column in arg_family_list])
    header.extend([f"Bacteria Defence Type '{column}' Count" for column in dg.columns])
    header.extend([f"Prophage Defence Type '{column}' Count" for column in dg.columns])
    header.extend([f"Bacteria Defence Subtype '{column}' Count" for column in defence_subtype_list])
    header.extend([f"Prophage Defence Subtype '{column}' Count" for column in defence_subtype_list])
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

with open('prophage_gene_impact_output.tsv', 'a') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    with open('results.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            df = pd.read_csv('results.tsv', sep='\t')
            biosampleid = df['BiosampleID']
            genome_id = biosampleid[row_num]
            prophage_arg_count = df[f"Prophage AMR Gene Count"][row_num]
            bacteria_arg_count = df[f"Bacteria ARG Count"][row_num]
            bacteria_defence_count = df[f"Bacteria Phage Defense System Count"][row_num]
            prophage_defence_count = df[f"Prophage Defense System Count"][row_num]
            bacteria_virulence_count = df[f"Bacteria Virulence Gene Count"][row_num]
            prophage_virulence_count = df[f"Prophage Virulence Gene Count"][row_num]
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

            bacteria_arg_family_counts = {}
            prophage_arg_family_counts = {}
            for column in arg_family_list:
                bacteria_arg_family_counts[column] = 0
                prophage_arg_family_counts[column] = 0
                bacteria_arg_family_counts[column] += df[f"Bacteria ARG Family '{column}' Count"][
                    row_num]
                column_name = f"Prophage ARG Family '{column}' Count"
                if column_name in df.columns:
                    prophage_arg_family_counts[column] += df[f"Prophage ARG Family '{column}' Count"][
                        row_num]


            bacteria_defence_subtype_counts = {}
            prophage_defence_subtype_counts = {}
            for column in defence_subtype_list:
                bacteria_defence_subtype_counts[column] = 0
                prophage_defence_subtype_counts[column] = 0
                bacteria_defence_subtype_counts[column] += df[f"Bacteria Phage Defense System Subtype '{column}' Count"][
                    row_num]
                column_name = f"Prophage Phage Defense System Subtype '{column}' Count"
                if column_name in df.columns:
                    prophage_defence_subtype_counts[column] += df[f"Prophage Phage Defense System Subtype '{column}' Count"][
                        row_num]

            row_num +=1


            output_row = [f'{genome_id}', f'{bacteria_arg_count}', f'{prophage_arg_count}', f'{bacteria_defence_count}', f'{prophage_defence_count}', f'{bacteria_virulence_count}', f'{prophage_virulence_count}']
            output_row.extend([f'{bacteria_arg_counts[column]}' for column in arg_columns])
            output_row.extend([f'{prophage_arg_counts[column]}' for column in arg_columns])
            output_row.extend([f'{bacteria_arg_family_counts[column]}' for column in arg_family_list])
            output_row.extend([f'{prophage_arg_family_counts[column]}' for column in arg_family_list])
            output_row.extend([f'{bacteria_defence_counts[column]}' for column in defence_columns])
            output_row.extend([f'{prophage_defence_counts[column]}' for column in defence_columns])
            output_row.extend([f'{bacteria_defence_subtype_counts[column]}' for column in defence_subtype_list])
            output_row.extend([f'{prophage_defence_subtype_counts[column]}' for column in defence_subtype_list])
            output_row.extend([f'{bacteria_virulence_counts[column]}' for column in virulence_list])
            output_row.extend([f'{prophage_virulence_counts[column]}' for column in virulence_list])
            tsv_writer.writerow(output_row)


# Read the output TSV file containing the summary
output_file = 'prophage_gene_impact_output.tsv'
output_df = pd.read_csv(output_file, sep='\t')

# Make output file to save results in
with open('prophage_gene_type_copies.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    header = ["Gene Name", "No Copies", "No Copies in Prophage", "Only Copies in Prophage", "1 Additional Copy",
              "2 Additional Copies", "3 Additional Copies", "4 Additional Copies", "5+ Additional Copies"]
    tsv_writer.writerow(header)

# Initialize a dictionary to store the counts
value_counts = {}

# Iterate over each column in the DataFrame
for column in output_df.columns[1:]:  # Exclude the 'BiosampleID' column
    if "Bacteria" in column:
        parts = column.split()
        # Take everything except the first word and join back to get name of column
        column_name = ' '.join(parts[1:])

        # Find matching prophage column
        prophage_col = None
        for c in output_df.columns[1:]:
            if c.startswith(f"Prophage {column_name}"):
                prophage_col = c
                break

        # Initialize count for different values
        no_copies = 0
        bacteria_only = 0
        prophage_only = 0
        one_additional = 0
        two_additional = 0
        three_additional = 0
        four_additional = 0
        five_plus_additional = 0
        bacteria_copies = 0
        prophage_copies = 0

        # Iterate over rows by index
        for idx in output_df.index:
            bacteria_copies = int(output_df.at[idx, column])
            prophage_copies = int(output_df.at[idx, prophage_col])

            if bacteria_copies == 0:
                no_copies += 1
            else:
                if prophage_copies == 0:
                    bacteria_only += 1
                elif bacteria_copies - prophage_copies == 0:
                    prophage_only += 1
                elif bacteria_copies - prophage_copies == 1:
                    one_additional += 1
                elif bacteria_copies - prophage_copies == 2:
                    two_additional += 1
                elif bacteria_copies - prophage_copies == 3:
                    three_additional += 1
                elif bacteria_copies - prophage_copies == 4:
                    four_additional += 1
                elif bacteria_copies - prophage_copies > 4:
                    five_plus_additional += 1

        with open('prophage_gene_type_copies.tsv', 'a') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            output_row = [f'{column_name}', f'{no_copies}', f'{bacteria_only}', f'{prophage_only}', f'{one_additional}',
                          f'{two_additional}', f'{three_additional}', f'{four_additional}', f'{five_plus_additional}']
            tsv_writer.writerow(output_row)
