import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import gzip
import socket
from Bio.SeqUtils import gc_fraction
import argparse
from argparse import ArgumentParser

parser = argparse.ArgumentParser()

parser.add_argument("cores", help="Number of Cores to use", type=int)
args = parser.parse_args()

# Write headings for file results for each genome will be recorded in
if not os.path.isfile('results.tsv'):
    with open('results.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ['BiosampleID', 'Genome Length', 'Genome Defense System Count', 'Genome AMR Count',
             'Genome Virulence Gene Count'])
# Write headings for file the virulence gene results will be recorded in for each genome
if not os.path.isfile('genome_virulence_genes.tsv'):
    with open('genome_virulence_genes.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Virulence Gene ID', 'Genome ID', 'Sequence', 'Start', 'End', 'Gene',
                             'Coverage', '% Coverage', '% Identity', 'Database', 'Accession', 'Gene Product'])
# Write headings for file on the AMR genes present in each Enterococcus genomes
if not os.path.isfile('genome_amr_genes.tsv'):
    with open('genome_amr_genes.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['AMR Gene ID', 'ORF ID', 'Source Sequence', 'Cut-Off', '% Identity', 'Drug Class',
                             'Resistance Mechanism', 'AMR Gene Family', '% Length of Reference Genome', 'Antibiotic',
                             'Predicted DNA', 'Predicted Protein'])
# Write headings for file on the phage defense systems in each Enterococcus genome
if not os.path.isfile('genome_phage_defenses.tsv'):
    with open('genome_phage_defenses.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['ID', 'Genome ID', 'System ID', 'Defense Type', 'Defense Subtype', 'System Start',
                             'System End', 'Proteins in System', 'Genes in System', 'Profiles in System'])

genome_defense_systems = 0
genome_amr_genes = 0
genome_virulence_genes = 0

for filename in os.listdir():
    if filename.startswith('SAM') and filename.endswith('.fna') and not filename.endswith('_rgi.fna'):
        genome_defense_systems = 0
        genome_amr_genes = 0
        genome_virulence_genes = 0
        id = filename.split('.')[0]
        if not os.path.isfile(f'{id}_rgi.txt'):
            # Get genome length
            genome_length = sum(len(r) for r in SeqIO.parse(f"{id}.fna", "fasta"))

            # Run defense-finder on enterococcus genome
            subprocess.run(f"defense-finder run {id}.fna", shell=True)

            with open(f'{id}_defense_finder_systems.tsv', 'r') as file:
                reader = csv.DictReader(file, delimiter='\t')
                row_num = 1

                # Record any defense systems found in the file
                for row in reader:
                    system_id = row['sys_id']
                    defense_type = row['type']
                    defense_subtype = row['subtype']
                    system_start = row['sys_beg']
                    system_end = row['sys_end']
                    proteins_in_system = row['protein_in_syst']
                    genes_in_system = row['genes_count']
                    profiles_in_system = row['name_of_profiles_in_sys']

                    with open('genome_phage_defenses.tsv', 'a') as out_file:
                        tsv_writer = csv.writer(out_file, delimiter='\t')
                        tsv_writer.writerow([f'{id}_sys{row_num}', f'{id}', f'{system_id}', f'{defense_type}',
                                             f'{defense_subtype}', f'{system_start}', f'{system_end}',
                                             f'{proteins_in_system}', f'{genes_in_system}', f'{profiles_in_system}'])
                    row_num +=1
                    genome_defense_systems +=1
            subprocess.run(f'rm *.prt*', shell=True)

            # Run rgi on enterococcus genome
            subprocess.run(f"rgi main -i {id}.fna -o {id}_rgi -n {args.cores} --include_loose --clean", shell=True)

            with open(f'{id}_rgi.txt', 'r') as file:
                reader = csv.DictReader(file, delimiter='\t')
                row_num = 1

                # Record any amr genes found in the genome_amr_genes.tsv file
                for row in reader:
                    amr_gene_id = f'{id}_amr{row_num}'
                    orf_id = row['ORF_ID']
                    source_sequence = row['Contig']
                    cut_off = row['Cut_Off']
                    percentage_coverage = row['Best_Identities']
                    drug_class = row['Drug Class']
                    resistance = row['Resistance Mechanism']
                    amr_gene_family = row['AMR Gene Family']
                    percentage_length = row['Percentage Length of Reference Sequence']
                    antibiotic = row['Antibiotic']
                    predicted_DNA = row['Predicted_DNA']
                    predicted_protein = row['Predicted_Protein']

                    if float(percentage_length) >= 75 and float(percentage_coverage) >= 90:
                        with open('genome_amr_genes.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{amr_gene_id}', f'{orf_id}', f'{source_sequence}', f'{cut_off}',
                                                 f'{percentage_coverage}' f'{drug_class}', f'{resistance}',
                                                 f'{amr_gene_family}', f'{percentage_length}', f'{antibiotic}',
                                                 f'{predicted_DNA}', f'{predicted_protein}'])
                        row_num += 1
                        genome_amr_genes += 1
                    elif cut_off != "Loose":
                        with open('genome_amr_genes.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{amr_gene_id}', f'{orf_id}', f'{source_sequence}', f'{cut_off}',
                                                 f'{percentage_coverage}' f'{drug_class}', f'{resistance}',
                                                 f'{amr_gene_family}', f'{percentage_length}', f'{antibiotic}',
                                                 f'{predicted_DNA}', f'{predicted_protein}'])
                        row_num += 1
                        genome_amr_genes += 1
                    else:
                        row_num += 1
                        # Run genome through abricate to identify virulence genes
            subprocess.run(f"abricate --db vfdb --minid 90 --mincov 75 {id}.fna > {id}_output.tsv",
                           shell=True)
            if os.path.isfile(f'{id}_output.tsv'):
                with open(f'{id}_output.tsv', 'r') as file:
                    reader = csv.DictReader(file, delimiter='\t')
                    row_num = 1

                    # Record any virulence genes found in the file
                    for row in reader:
                        genome_id = row['#FILE']
                        sequence = row['SEQUENCE']
                        sequence_start = row['START']
                        sequence_end = row['END']
                        gene_name = row['GENE']
                        coverage = row['COVERAGE']
                        percent_coverage = row['%COVERAGE']
                        percent_identity = row['%IDENTITY']
                        database = row['DATABASE']
                        vg_accession = row['ACCESSION']
                        gene_product = row['PRODUCT']


                        with open('genome_virulence_genes.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow(
                                [f'{id}__{row_num}', f'{genome_id}', f'{sequence}',
                                 f'{sequence_start}', f'{sequence_end}', f'{gene_name}', f'{coverage}',
                                 f'{percent_coverage}', f'{percent_identity}', f'{database}', f'{vg_accession}',
                                 f'{gene_product}'])
                        row_num += 1
                        genome_virulence_genes += 1

            # Add all information to results.tsv, for analysis, only want first two words in Organism
            with open('results.tsv', 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow([f'{id}', f'{genome_length}', f'{genome_defense_systems}', f'{genome_amr_genes}',
                                     f'{genome_virulence_genes}'])
