import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import gzip

if not os.path.isdir('Output'):
    os.mkdir('Output')

df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid']

# Write headings for file results for each genome will be recorded in
# Table containing list of all the results for each prophage
if not os.path.isfile('Output/results.tsv'):
    with open('Output/results.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ['Biosampleid',
             'Number of Prophages Containing Phage Defences', 'Total Number of Prophage Phage Defences',
             'Number of Prophages Containing Anti-defences', 'Total Number of Prophage Anti-defences'])
if not os.path.isfile('Output/prophage_results.tsv'):
    with open('Output/prophage_results.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ['Prophage Name', 'Biosampleid', 'Number of Phage Defences', 'Number of Antidefences'])
# Write headings for file the rgi-card results will be recorded in for each genome
if not os.path.isfile('Output/antidefense.tsv'):
    with open('Output/antidefense.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['ID', 'Genome ID', 'System ID', 'Defense Type', 'Defense Subtype', 'System Start',
                             'System End', 'Proteins in System', 'Genes in System', 'Profiles in System'])
if not os.path.isfile('Output/phagedefense.tsv'):
    with open('Output/phagedefense.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['ID', 'Genome ID', 'System ID', 'Defense Type', 'Defense Subtype', 'System Start',
                             'System End', 'Proteins in System', 'Genes in System', 'Profiles in System'])

for id in biosampleid:
    if os.path.isfile(f'{id}.tar.gz'):
        subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)  # Unzip the file for the biosampleID

        os.chdir(f'{id}/phageboost_results')

        prophage_num = 0
        antidefenses = 0
        antidefense_prophages = 0
        phagedefenses = 0
        phagedefense_prophages = 0

        file_num = 1

        for filename in os.listdir():
            if filename.startswith(f'{id}_p'):
                prophage_num += 1  # Get number of prophages, starts at 0, adds 1 each cycle for each directory
                os.chdir(f'{id}_p{file_num}')

                # Run defense-finder on enterococcus genome
                subprocess.run(f"defense-finder run {id}_p{file_num}.fna -o phagedefenses", shell=True)
                subprocess.run(f'rm phagedefenses/*.prt*', shell=True)
                with open(f'phagedefenses/{id}_p{file_num}_defense_finder_systems.tsv', 'r') as file:
                    reader = csv.DictReader(file, delimiter='\t')
                    row_num = 1
                    prophage_defenses = 0

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

                        os.chdir('../../../Output')
                        with open('phagedefense.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{id}_sys{row_num}', f'{id}', f'{system_id}', f'{defense_type}',
                                                 f'{defense_subtype}', f'{system_start}', f'{system_end}',
                                                 f'{proteins_in_system}', f'{genes_in_system}', f'{profiles_in_system}'])
                        os.chdir('..')
                        row_num += 1
                        prophage_phage_defenses += 1
                        phagedefenses +=1
                        os.chdir(f'{id}/phageboost_results/{id}_p{file_num}')

                # Run defense-finder on enterococcus genome
                subprocess.run(f"defense-finder run {id}_p{file_num}.fna -A -o antidefenses", shell=True)
                subprocess.run(f'rm antidefenses/*.prt*', shell=True)
                with open(f'antidefenses/{id}_p{file_num}_defense_finder_systems.tsv', 'r') as file:
                    reader = csv.DictReader(file, delimiter='\t')
                    row_num = 1
                    prophage_antidefenses = 0

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

                        os.chdir('../../../Output')
                        with open('antidefense.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{id}_sys{row_num}', f'{id}', f'{system_id}', f'{defense_type}',
                                                 f'{defense_subtype}', f'{system_start}', f'{system_end}',
                                                 f'{proteins_in_system}', f'{genes_in_system}', f'{profiles_in_system}'])
                        os.chdir('..')
                        row_num += 1
                        prophage_antidefenses += 1
                        antidefenses +=1
                        os.chdir(f'{id}/phageboost_results/{id}_p{file_num}')

                os.chdir('../../../Output')
                with open('prophage_results.tsv', 'a') as out_file:
                    tsv_writer = csv.writer(out_file, delimiter='\t')
                    tsv_writer.writerow(
                        [f'{id}_p{file_num}', f'{prophage_defenses}', f'{prophage_antidefenses}'])

                os.chdir(f'../{id}/phageboost_results')
                file_num += 1
                if prophage_defenses > 0:
                    phagedefense_prophages += 1
                if prophage_antidefenses > 0:
                    antidefense_prophages += 1

        os.chdir('../../Output')

        # Add all information to results.tsv, for analysis, only want first two words in Organism
        with open('results.tsv', 'a') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            tsv_writer.writerow([f'{id}', f'{phagedefense_prophages}', f'{phagedefenses}',
                                 f'{antidefense_prophages}', f'{antidefenses}'])
        os.chdir('..')
        subprocess.run(f"tar -zcf {id}.tar.gz {id}", shell=True)
        shutil.rmtree(id)
        print(f"{id} Done")

os.chdir('Output')