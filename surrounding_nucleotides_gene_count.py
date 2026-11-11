import os
import pandas as pd
import csv

before_args = pd.read_csv('before_args.tab', sep='\t')
before_virs = pd.read_csv('before_virulence_genes.tab', sep='\t')
before_defs = pd.read_csv('prophage_nucleotides_before_defense_finder_systems.tsv', sep='\t')
before_anti_defs = pd.read_csv('prophage_nucleotides_before_anti_defense_finder_systems.tsv', sep='\t')
after_args = pd.read_csv('after_args.tab', sep='\t')
after_virs = pd.read_csv('after_virulence_genes.tab', sep='\t')
after_defs = pd.read_csv('prophage_nucleotides_after_defense_finder_systems.tsv', sep='\t')
after_anti_defs = pd.read_csv('prophage_nucleotides_after_anti_defense_finder_systems.tsv', sep='\t')

before_virulence_list = before_virs['GENE'].unique()
after_virulence_list = after_virs['GENE'].unique()
before_args_list = before_args['GENE'].unique()
after_args_list = after_args['GENE'].unique()
before_defs_list = before_defs['subtype'].unique()
before_anti_defs_list = before_anti_defs['subtype'].unique()
after_defs_list = after_defs['subtype'].unique()
after_anti_defs_list = after_anti_defs['subtype'].unique()

with open('surrounding_nucleotides_gene_count.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    header = ["Motif", "Before Virulence Gene Count", "Before ARG Count", "Before Defence System Count",
              "Before Anti-Defence System Count", "After Virulence Gene Count", "After ARG Count",
              "After Defence System Count", "After Anti-Defence System Count"]
    header.extend([f"Before Virulence Gene '{column}' Count" for column in before_virulence_list])
    header.extend([f"After Virulence Gene '{column}' Count" for column in after_virulence_list])
    header.extend([f"Before ARG '{column}' Count" for column in before_args_list])
    header.extend([f"After ARG '{column}' Count" for column in after_args_list])
    header.extend([f"Before Defence Subtype '{column}' Count" for column in before_defs_list])
    header.extend([f"After Defence Subtype '{column}' Count" for column in after_defs_list])
    header.extend([f"Before Anti-Defence Subtype '{column}' Count" for column in before_anti_defs_list])
    header.extend([f"After Anti-Defence Subtype '{column}' Count" for column in after_anti_defs_list])
    tsv_writer.writerow(header)


row_num = 0
before_motif_number = 1


while before_motif_number < 22855:
    before_arg_count = 0
    after_arg_count = 0
    before_defence_count = 0
    after_defence_count = 0
    before_virulence_count = 0
    after_virulence_count = 0
    before_anti_defence_count = 0
    after_anti_defence_count = 0
    before_virulence_counts = {}
    for column in before_virulence_list:
        before_virulence_counts[column] = 0
    after_virulence_counts = {}
    for column in after_virulence_list:
        after_virulence_counts[column] = 0
    before_arg_counts = {}
    for column in before_args_list:
        before_arg_counts[column] = 0
    after_arg_counts = {}
    for column in after_args_list:
        after_arg_counts[column] = 0
    before_def_counts = {}
    for column in before_defs_list:
        before_def_counts[column] = 0
    after_def_counts = {}
    for column in after_defs_list:
        after_def_counts[column] = 0
    before_anti_def_counts = {}
    for column in before_anti_defs_list:
        before_anti_def_counts[column] = 0
    after_anti_def_counts = {}
    for column in after_anti_defs_list:
        after_anti_def_counts[column] = 0
    with open('before_args.tab', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f"motif_{before_motif_number}_before" in row["SEQUENCE"]:
                before_arg_count += 1
                gene = row["GENE"]
                before_arg_counts[gene] += 1
    with open('before_virulence_genes.tab', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f"motif_{before_motif_number}_before" in row["SEQUENCE"]:
                before_virulence_count += 1
                gene = row["GENE"]
                before_virulence_counts[gene] += 1
    with open('prophage_nucleotides_before_defense_finder_systems.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f"motif_{before_motif_number}_before" in row['sys_id']:
                before_defence_count += 1
                subtype = row["subtype"]
                before_def_counts[subtype] += 1
    with open('prophage_nucleotides_before_anti_defense_finder_systems.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f"motif_{before_motif_number}_before" in row['sys_id']:
                before_anti_defence_count += 1
                subtype = row["subtype"]
                before_anti_def_counts[subtype] += 1
    with open('surrounding_nucleotides_gene_count.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        output_row = [f"before_motif_{before_motif_number}", before_virulence_count, before_arg_count,
                      before_defence_count, before_anti_defence_count, after_virulence_count, after_arg_count,
                      after_defence_count,after_anti_defence_count]
        output_row.extend(f"{before_virulence_counts[column]}" for column in before_virulence_list)
        output_row.extend(f"{after_virulence_counts[column]}" for column in after_virulence_list)
        output_row.extend(f"{before_arg_counts[column]}" for column in before_args_list)
        output_row.extend(f"{after_arg_counts[column]}" for column in after_args_list)
        output_row.extend(f"{before_def_counts[column]}" for column in before_defs_list)
        output_row.extend(f"{after_def_counts[column]}" for column in after_defs_list)
        output_row.extend(f"{before_anti_def_counts[column]}" for column in before_anti_defs_list)
        output_row.extend(f"{after_anti_def_counts[column]}" for column in after_anti_defs_list)
        tsv_writer.writerow(output_row)
    before_motif_number += 1

after_motif_number = 1

while after_motif_number < 28898:
    before_arg_count = 0
    after_arg_count = 0
    before_defence_count = 0
    after_defence_count = 0
    before_virulence_count = 0
    after_virulence_count = 0
    before_anti_defence_count = 0
    after_anti_defence_count = 0
    before_virulence_counts = {}
    for column in before_virulence_list:
        before_virulence_counts[column] = 0
    after_virulence_counts = {}
    for column in after_virulence_list:
        after_virulence_counts[column] = 0
    before_arg_counts = {}
    for column in before_args_list:
        before_arg_counts[column] = 0
    after_arg_counts = {}
    for column in after_args_list:
        after_arg_counts[column] = 0
    before_def_counts = {}
    for column in before_defs_list:
        before_def_counts[column] = 0
    after_def_counts = {}
    for column in after_defs_list:
        after_def_counts[column] = 0
    before_anti_def_counts = {}
    for column in before_anti_defs_list:
        before_anti_def_counts[column] = 0
    after_anti_def_counts = {}
    for column in after_anti_defs_list:
        after_anti_def_counts[column] = 0
    with open('after_args.tab', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f"motif_{after_motif_number}_after" in row["SEQUENCE"]:
                after_arg_count += 1
                gene = row["GENE"]
                after_arg_counts[gene] += 1
    with open('after_virulence_genes.tab', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f"motif_{after_motif_number}_after" in row["SEQUENCE"]:
                after_virulence_count += 1
                gene = row["GENE"]
                after_virulence_counts[gene] += 1
    with open('prophage_nucleotides_after_defense_finder_systems.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f"motif_{after_motif_number}_after" in row['sys_id']:
                after_defence_count += 1
                subtype = row["subtype"]
                after_def_counts[subtype] += 1
    with open('prophage_nucleotides_after_anti_defense_finder_systems.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f"motif_{after_motif_number}_after" in row['sys_id']:
                after_anti_defence_count += 1
                subtype = row["subtype"]
                after_anti_def_counts[subtype] += 1
    with open('surrounding_nucleotides_gene_count.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        output_row = [f"after_motif_{after_motif_number}", before_virulence_count, before_arg_count,
                      before_defence_count, before_anti_defence_count, after_virulence_count, after_arg_count,
                      after_defence_count,after_anti_defence_count]
        output_row.extend(f"{before_virulence_counts[column]}" for column in before_virulence_list)
        output_row.extend(f"{after_virulence_counts[column]}" for column in after_virulence_list)
        output_row.extend(f"{before_arg_counts[column]}" for column in before_args_list)
        output_row.extend(f"{after_arg_counts[column]}" for column in after_args_list)
        output_row.extend(f"{before_def_counts[column]}" for column in before_defs_list)
        output_row.extend(f"{after_def_counts[column]}" for column in after_defs_list)
        output_row.extend(f"{before_anti_def_counts[column]}" for column in before_anti_defs_list)
        output_row.extend(f"{after_anti_def_counts[column]}" for column in after_anti_defs_list)
        tsv_writer.writerow(output_row)
    after_motif_number += 1
