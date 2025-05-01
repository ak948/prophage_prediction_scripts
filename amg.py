import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import gzip

working_directory = os.getcwd()
if not os.path.isdir('amgOutput'):
    os.mkdir('amgOutput')

df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid']

# Write headings for file results for each genome will be recorded in
# Table containing list of all the results for each prophage
if not os.path.isfile('amgOutput/results.tsv'):
    with open('amgOutput/results.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ['Biosampleid', 'Number of Prophages Containing AMGs', 'Total Number of Prophage AMGs',
             'Total Number of Prophage KEGG Pathways'])
if not os.path.isfile('amgOutput/prophage_results.tsv'):
    with open('amgOutput/prophage_results.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ['Prophage Name', 'Biosampleid', 'Number of AMGs', 'Number of KEGG Pathways'])
if not os.path.isfile('amgOutput/amg_pathways.tsv'):
    with open('amgOutput/amg_pathways.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ['ID', 'Prophage ID', 'BiosampleID', 'KEGG Entry', 'Metabolism', 'Pathways', 'Total AMGs',
             'Present AMG KOs'])
# Write headings for file the rgi-card results will be recorded in for each genome
if not os.path.isfile('amgOutput/amg_individuals.tsv'):
    with open('amgOutput/amg_individuals.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['ID', 'Prophage ID', 'BiosampleID', 'Protein', 'AMG KO', 'AMG KO Name', 'Pfam',
                             'Pfam Name'])

prophage_amgs = 0
amg_prophages = 0
kegg_pathways = 0

for id in biosampleid:
    if os.path.isfile(f'{id}.tar.gz'):
        subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)  # Unzip the file for the biosampleID

        os.chdir(f'{id}/phageboost_results')

        prophage_num = 0
        prophage_amgs = 0
        amg_prophages = 0
        kegg_pathways = 0

        file_num = 1

        for filename in os.listdir():
            if filename.startswith(f'{id}_p'):
                prophage_num += 1  # Get number of prophages, starts at 0, adds 1 each cycle for each directory
                os.chdir(f'{id}_p{file_num}')

                # Run defense-finder on enterococcus genome
                subprocess.run(f"python3 /shared/team/VIBRANT/VIBRANT_run.py -i {id}_p{file_num}.fna -d /shared/team/conda/andym.millardlab/akinsphageboost/share/vibrant-1.2.1/db/databases/ -t 40 -no_plot", shell=True)
                os.chdir(f'VIBRANT_{id}_p{file_num}/VIBRANT_results_{id}_p{file_num}')
                with open(f'VIBRANT_AMG_individuals_{id}_p{file_num}.tsv', 'r') as file:
                    reader = csv.DictReader(file, delimiter='\t')
                    row_num = 1
                    amgs = 1

                    # Record any defense systems found in the file
                    for row in reader:
                        protein = row['protein']
                        amg_ko = row['AMG KO']
                        amg_ko_name = row['AMG KO name']
                        pfam = row['Pfam']
                        pfam_name = row['Pfam name']

                        os.chdir(f'{working_directory}/amgOutput')
                        with open('amg_individuals.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{id}_p{file_num}_amg{row_num}', f'{id}_p{file_num}', f'{id}',
                                                 f'{protein}', f'{amg_ko}', f'{amg_ko_name}', f'{pfam}',
                                                 f'{pfam_name}'])
                        row_num +=1
                        amgs +=1
                        os.chdir(f'{working_directory}/{id}/phageboost_results/{id}_p{file_num}')

                os.chdir(f'{working_directory}/{id}/phageboost_results/{id}_p{file_num}/VIBRANT_{id}_p{file_num}/VIBRANT_results_{id}_p{file_num}')
                with open(f'VIBRANT_AMG_pathways_{id}_p{file_num}.tsv', 'r') as file:
                    reader = csv.DictReader(file, delimiter='\t')
                    row_num = 1
                    keggs = 1

                    # Record any defense systems found in the file
                    for row in reader:
                        kegg_entry = row['KEGG Entry']
                        metabolism = row['Metabolism']
                        pathway = row['Pathway']
                        total_amgs = row['Total AMGs']
                        present_amg_kos = row['Present AMG KOs']

                        os.chdir(f'{working_directory}/amgOutput')
                        with open('amg_pathways.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{id}_p{file_num}_kegg{row_num}', f'{id}_p{file_num}', f'{id}',
                                                 f'{kegg_entry}', f'{metabolism}', f'{pathway}', f'{total_amgs}',
                                                 f'{present_amg_kos}'])
                        os.chdir('..')
                        row_num += 1
                        keggs +=1
                        os.chdir(f'{id}/phageboost_results/{id}_p{file_num}')

                os.chdir(f'{working_directory}/{id}/phageboost_results/{id}_p{file_num}')
                subprocess.run(f'rm -r VIBRANT*', shell=True)
                os.chdir('../../../amgOutput')
                with open('prophage_results.tsv', 'a') as out_file:
                    tsv_writer = csv.writer(out_file, delimiter='\t')
                    tsv_writer.writerow(
                        [f'{id}_p{file_num}', f'{amgs}',f'{keggs}'])

                os.chdir(f'../{id}/phageboost_results')
                file_num += 1
                if amgs > 0:
                    amg_prophages += 1
                    prophage_amgs += amgs
                    kegg_pathways += keggs

        os.chdir('../../amgOutput')

        # Add all information to results.tsv, for analysis, only want first two words in Organism
        with open('results.tsv', 'a') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            tsv_writer.writerow([f'{id}', f'{amg_prophages}', f'{prophage_amgs}', f'{kegg_pathways}'])
        os.chdir('..')
        shutil.rmtree(id)
        print(f"{id} Done")

os.chdir('amgOutput')