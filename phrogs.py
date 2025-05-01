import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import shutil
import numpy as np
import gzip
import sys

csv.field_size_limit(sys.maxsize)

with open('PhrogOutput/phrog_counts.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Prophage ID', 'PHROGs Genes Count', 'Non-PHROGs Genes Count', '% of PHROGs Genes',
                         'PHROGs Modification Genes Count', 'PHROGs Resistance Genes Count',
                         'PHROGs Virulence Genes Count', 'PHROGs Toxin Genes Count', 'PHROGs Immune Genes Count',
                         'PHROGs Defense Genes Count', 'Modification and Toxin Genes Count',
                         'Modification and Defense Genes Count', 'Resistance and Toxin Genes Count',
                         'Virulence and Toxin Genes Count'])

with open('PhrogOutput/phrog_genes.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Gene ID', 'Prophage ID', 'PHROGs Gene Category', 'PHROGs Gene ID'])

df = pd.read_csv('results.tsv', sep='\t')
biosampleid = df['BiosampleID']

working_directory = os.getcwd()

for id in biosampleid:
    os.chdir(working_directory)
    if os.path.isfile(f'{id}.tar.gz'):
        subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)  # Unzip the file for the biosampleID
        os.chdir(f'{id}/phageboost_results')

        file_num = 1

        for filename in os.listdir():
            if filename.startswith(f'{id}_p'):
                os.chdir(f'{id}_p{file_num}')
                result = subprocess.run(f"cat {id}_p{file_num}.gbk | grep 'all_phrogs:' | wc -l", shell = True,stdout=subprocess.PIPE, text=True)
                phrogs = int(result.stdout)
                genes = 0
                for seq_record in SeqIO.parse(f"{id}_p{file_num}.gbk", "genbank"):  # Get number of CDS in the prophage
                    genes = len(seq_record.features)
                    genes -= 1
                non_phrogs = int(genes) - int(phrogs)
                phrogs_percent = int(phrogs) / int(genes) * 100

                # List of modification genes
                dt = pd.read_csv(f'{working_directory}/modification_genes.txt', sep='\t')
                modification_gene_list = dt['phrog']
                modification_gene_count = 0
                modification_and_toxin_gene_count = 0
                modification_and_defense_gene_count = 0

                # List of resistance genes
                dg = pd.read_csv(f'{working_directory}/resistance_genes.txt', sep='\t')
                resistance_gene_list = dg['phrog']
                resistance_gene_count = 0
                resistance_and_toxin_gene_count = 0

                # List of virulence genes
                dh = pd.read_csv(f'{working_directory}/virulence_genes.txt', sep='\t')
                virulence_gene_list = dh['phrog']
                virulence_gene_count = 0
                virulence_and_toxin_gene_count = 0

                # List of toxin genes
                di = pd.read_csv(f'{working_directory}/toxin_genes.txt', sep='\t')
                toxin_gene_list = di['phrog']
                toxin_gene_count = 0

                # List of immune genes
                dj = pd.read_csv(f'{working_directory}/immune_genes.txt', sep='\t')
                immune_gene_list = dj['phrog']
                immune_gene_count = 0

                # List of defense genes
                dk = pd.read_csv(f'{working_directory}/defense_genes.txt', sep='\t')
                defense_gene_list = dk['phrog']
                defense_gene_count = 0

                # Run the command and get the output
                output = subprocess.check_output(f"cat {id}_p{file_num}.gbk | grep 'all_phrogs:' || true", shell=True)

                if output:
                    # Split the output into lines
                    lines = output.decode().split('\n')

                    gene_num = 0
                    # Iterate over each line
                    for line in lines:
                        if 'all_phrogs:' in line:
                            phrog_number = line.split(':')[2].split('"')[0]

                            # Check if phrog_number matches the list of saved phrog_numbers
                            if phrog_number in modification_gene_list.values and phrog_number in toxin_gene_list.values:
                                gene_category = 'Modification and Toxin Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num += 1
                                modification_and_toxin_gene_count += 1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in modification_gene_list.values and phrog_number in defense_gene_list.values:
                                gene_category = 'Modification and Defense Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num += 1
                                modification_and_defense_gene_count += 1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in resistance_gene_list.values and phrog_number in toxin_gene_list.values:
                                gene_category = 'Resistance and Toxin Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num += 1
                                resistance_and_toxin_gene_count += 1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in virulence_gene_list.values and phrog_number in toxin_gene_list.values:
                                gene_category = 'Virulence and Toxin Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num += 1
                                virulence_and_toxin_gene_count += 1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in modification_gene_list.values:
                                gene_category = 'Modification Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num +=1
                                modification_gene_count +=1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in resistance_gene_list.values:
                                gene_category = 'Resistance Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num +=1
                                resistance_gene_count +=1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in virulence_gene_list.values:
                                gene_category = 'Virulence Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num +=1
                                virulence_gene_count +=1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in toxin_gene_list.values:
                                gene_category = 'Toxin Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num +=1
                                toxin_gene_count +=1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in immune_gene_list.values:
                                gene_category = 'Immune Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num +=1
                                immune_gene_count +=1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])
                            elif phrog_number in defense_gene_list.values:
                                gene_category = 'Defense Gene'
                                gene_id = f'{id}_p{file_num}_gene{gene_num}'
                                gene_num +=1
                                defense_gene_count +=1
                                with open(f'{working_directory}/PhrogOutput/phrog_genes.tsv', 'a') as out_file:
                                    tsv_writer = csv.writer(out_file, delimiter='\t')
                                    tsv_writer.writerow([f'{gene_id}', f'{id}', f'{gene_category}', f'{phrog_number}'])


                with open(f'{working_directory}/PhrogOutput/phrog_counts.tsv', 'a') as out_file:
                    tsv_writer = csv.writer(out_file, delimiter='\t')
                    tsv_writer.writerow([f'{id}_p{file_num}', f'{phrogs}', f'{non_phrogs}', f'{phrogs_percent}',
                                         f'{modification_gene_count}', f'{resistance_gene_count}',
                                         f'{virulence_gene_count}', f'{toxin_gene_count}', f'{immune_gene_count}',
                                         f'{defense_gene_count}', f'{modification_and_defense_gene_count}',
                                         f'{resistance_and_toxin_gene_count}', f'{resistance_and_toxin_gene_count}',
                                         f'{virulence_and_toxin_gene_count}'])
                os.chdir('..')
                file_num +=1
        os.chdir(working_directory)
        subprocess.run(f'rm -r {id}', shell=True)



