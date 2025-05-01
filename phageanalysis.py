import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import pylab
import matplotlib.pyplot as plt
from pylab import MaxNLocator
import gzip
from Bio import Entrez
import urllib
import time
from datetime import datetime
import socket
from Bio.SeqUtils import gc_fraction
import argparse
from argparse import ArgumentParser

parser = argparse.ArgumentParser()

parser.add_argument("cores", help="Number of Cores to use", type=int)
parser.add_argument("--start", help="Number of BiosampleID to start on (For splitting into chunks)", type=int)
parser.add_argument("--end", help="Number of BiosampleID to end on (For splitting into chunks)", type=int)
parser.add_argument("--skip_figures", help="Should figures be auto-produced at the end of the script", action="store_true")
args = parser.parse_args()

if not os.path.isdir('Output'):
    os.mkdir('Output')
# Get list of biosampleIDs
if type(args.start) is int:
    df = pd.read_csv('genomeinfo.tsv', sep='\t')
    biosampleid = df['Biosampleid'].loc[args.start : args.end]
else:
    df = pd.read_csv('genomeinfo.tsv', sep='\t')
    biosampleid = df['Biosampleid']

# Write headings for file results for each genome will be recorded in
if not os.path.isfile('Output/results.tsv'):
    with open('Output/results.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ['Biosampleid', 'Genbank Info Species', 'Collection Date', 'Collection Year', 'Collection Location',
             'Collection Country', 'Isolation Source', 'Number of Prophages', 'Average Prophage Length',
             'Minimum Prophage Length', 'Maximum Prophage Length',
             'Average number of prophage coding sequences', 'Minimum number of prophage coding sequences',
             'Maximum number of prophage coding sequences', 'Average Prophage GC%', 'Minimum Prophage GC%',
             'Maximum Prophage GC%', 'Number of prophages containing virulence genes',
             'Total number of virulence genes in prophages', 'Number of prophages containing AMR Genes',
             'Total number of AMR Genes in prophages'])
# Table containing list of all the results for each prophage
if not os.path.isfile('Output/prophage_results.tsv'):y
    with open('Output/prophage_results.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ['Prophage Name', 'Biosampleid', 'Genome Length', 'Coding Sequences', 'Virulence Gene Count', 'AMR Count', 'GC%'])
# Write headings for file the rgi-card results will be recorded in for each genome
if not os.path.isfile('Output/amr_genes.tsv'):
    with open('Output/amr_genes.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['AMR Gene ID', 'ORF ID', 'Source Sequence', 'Cut-Off', '% Identity', 'Drug Class',
                             'Resistance Mechanism', 'AMR Gene Family', '% Length of Reference Genome', 'Nudged', 'Note',
                             'Antibiotic', 'Predicted DNA', 'Predicted Protein'])
# Write headings for file the virulence gene results will be recorded in for each genome
if not os.path.isfile('Output/virulence_genes.tsv'):
    with open('Output/virulence_genes.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Virulence Gene ID', 'Prophage ID', 'Sequence', 'Start', 'End', 'Gene',
                             'Coverage', '% Coverage', '% Identity', 'Database', 'Accession', 'Gene Product'])
if os.path.isfile('Output/all_prophage_proteins.faa'):
    subprocess.run("rm all_prophage_proteins.faa", shell=True)
if not os.path.isdir('Output/fnasequences'):
    subprocess.run("mkdir Output/fnasequences", shell=True)


total_lengths = []
total_genes = []
total_amr_genes = []
total_virulence_genes = []


for id in biosampleid:
    subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)  # Unzip the file for the biosampleID
    # Want to take organism, collection date, collection location, isolation source from genomeinfo.tsv
    with open('genomeinfo.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        # Iterate over each row in the file
        for row in reader:
            if row['Biosampleid'] == id:
                # Extract the desired data from different columns
                organism = row['Organism']
                # Splitting the organism string to only keep the species
                split_organism = organism.split()
                enterococcus_index = split_organism.index("Enterococcus")

                # Extracting the species
                species = split_organism[enterococcus_index + 1]

                # Joining the words before Enterococcus
                words_before = " ".join(split_organism[:enterococcus_index + 1])
                final = words_before + ' ' + species

                isolation_date = row['Collection Date']
                date = row['Collection Date']
                if "-" in date:
                    split_date = date.split("-")
                    year = split_date[0]
                elif "/" in date:
                    split_date = date.split("/")
                    try:
                        year = split_date[2]
                    except:
                        year = split_date[1]
                else:
                    year = row['Collection Date']
                isolation_source = row['Isolation Source']
                isolation_location = row['Collection Location']
                split_country = isolation_location.split(":")
                country = split_country[0]

    os.chdir(f'{id}/phageboost_results')

    prophage_num = 0
    amr_gene_prophages = 0
    amr_genes = 0
    virulence_gene_prophages = 0
    virulence_genes = 0
    length_list = []
    genes_list = []
    gc_list = []

    file_num = 1

    for filename in os.listdir():
        if filename.startswith(f'{id}_p'):
            prophage_num += 1  # Get number of prophages, starts at 0, adds 1 each cycle for each directory
            os.chdir(f'{id}_p{file_num}')
            for seq_record in SeqIO.parse(f"{id}_p{file_num}.fna", "fasta"):  # Get sequence length of the prophage
                length = len(seq_record)
                length_list.append(length)
                total_lengths.append(int(length))
            for seq_record in SeqIO.parse(f"{id}_p{file_num}.gbk", "genbank"):  # Get number of CDS in the prophage
                genes = len(seq_record.features)
                genes_list.append(genes)
                total_genes.append(int(genes))

            gc_values = sorted(
                100 * gc_fraction(rec.seq) for rec in SeqIO.parse(f"{id}_p{file_num}.fna", "fasta"))  # Gene GC%
            gc_list.append(gc_values)

            # Run genome through abricate to identify virulence genes
            subprocess.run(f"abricate --db vfdb --minid 90 --mincov 75 {id}_p{file_num}.fna > output.tsv",
                           shell=True)
            if os.path.isfile('output.tsv'):
                with open('output.tsv', 'r') as file:
                    reader = csv.DictReader(file, delimiter='\t')
                    row_num = 1
                    prophage_virulence_genes = 0

                    # Record any virulence genes found in the file
                    for row in reader:
                        prophage_id = row['#FILE']
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

                        os.chdir('../../../Output')
                        with open('virulence_genes.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{id}_p_{file_num}_{row_num}', f'{prophage_id}', f'{sequence}',
                                                 f'{sequence_start}', f'{sequence_end}', f'{gene_name}', f'{coverage}',
                                                 f'{percent_coverage}', f'{percent_identity}', f'{database}', f'{vg_accession}',
                                                 f'{gene_product}'])
                        os.chdir('..')
                        row_num += 1
                        virulence_genes += 1
                        prophage_virulence_genes += 1
                        os.chdir(f'{id}/phageboost_results/{id}_p{file_num}')
            else:
                os.chdir('..')
                os.chdir(f'{id}/phageboost_results/{id}_p{file_num}')
                prophage_virulence_genes = 0


            # Run genome through RGI-CARD to identify AMR genes
            subprocess.run(f"rgi main -i {id}_p{file_num}.fna -o {id}_p{file_num}_rgi -n {args.cores} --include_loose --clean",
                           shell=True)
            os.remove(f'{id}_p{file_num}_rgi.json')
            with open(f'{id}_p{file_num}_rgi.txt', 'r') as file:
                reader = csv.DictReader(file, delimiter='\t')
                row_num = 1
                prophage_amr_genes = 0

                # Record any amr genes found in the amr_genes.tsv file
                for row in reader:
                    amr_gene_id = f'{id}_p{file_num}_amr{row_num}'
                    orf_id = row['ORF_ID']
                    source_sequence = row['Contig']
                    cut_off = row['Cut_Off']
                    percentage_coverage = row['Best_Identities']
                    drug_class = row['Drug Class']
                    resistance = row['Resistance Mechanism']
                    amr_gene_family = row['AMR Gene Family']
                    percentage_length = row['Percentage Length of Reference Sequence']
                    nudged = row['Nudged']
                    note = row['Note']
                    antibiotic = row['Antibiotic']
                    predicted_DNA = row['Predicted_DNA']
                    predicted_protein = row['Predicted_Protein']

                    if float(percentage_length) >= 75 and float(percentage_coverage) >= 90:
                        os.chdir('../../../Output')
                        with open('amr_genes.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{amr_gene_id}', f'{orf_id}', f'{source_sequence}', f'{cut_off}',
                                                 f'{percentage_coverage}', f'{drug_class}', f'{resistance}',
                                                 f'{amr_gene_family}',
                                                 f'{percentage_length}', f'{nudged}', f'{note}', f'{antibiotic}',
                                                 f'{predicted_DNA}', f'{predicted_protein}'])
                        os.chdir('..')
                        row_num += 1
                        amr_genes += 1
                        prophage_amr_genes += 1
                        os.chdir(f'{id}/phageboost_results/{id}_p{file_num}')
                    elif cut_off != "Loose":
                        os.chdir('../../../Output')
                        with open('amr_genes.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow([f'{amr_gene_id}', f'{orf_id}', f'{source_sequence}', f'{cut_off}',
                                                 f'{percentage_coverage}', f'{drug_class}', f'{resistance}',
                                                 f'{amr_gene_family}',
                                                 f'{percentage_length}', f'{nudged}', f'{note}', f'{antibiotic}',
                                                 f'{predicted_DNA}', f'{predicted_protein}'])
                        os.chdir('..')
                        row_num += 1
                        amr_genes += 1
                        prophage_amr_genes += 1
                        os.chdir(f'{id}/phageboost_results/{id}_p{file_num}')
                    else:
                        row_num += 1

            subprocess.run(f'rm {id}_p{file_num}_rgi.txt', shell=True)
            subprocess.run(f'rm output.tsv', shell=True)
            # Make .faa file containing nucleotide sequences of all prophages
            subprocess.run(f"cat {id}_p{file_num}.fna >> ../../../Output/all_prophage_sequences.fna", shell=True)
            subprocess.run(f"cp {id}_p{file_num}.fna ../../../Output/fnasequences", shell=True)

            os.chdir('../../../Output')
            with open('prophage_results.tsv', 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow(
                    [f'{id}_p{file_num}', f'{id}', f'{length}', f'{genes}', f'{prophage_virulence_genes}',
                     f'{prophage_amr_genes}', statistics.mean(gc_values)])

            os.chdir(f'../{id}/phageboost_results')
            file_num += 1
            if prophage_amr_genes > 0:
                amr_gene_prophages += 1
            if prophage_virulence_genes > 0:
                virulence_gene_prophages += 1
        total_amr_genes.append(int(amr_genes))
        total_virulence_genes.append(int(virulence_genes))

    os.chdir('../../Output')
    if prophage_num != 0:
        # Convert lengths and genes into numpy arrays
        lengths = np.array(length_list)
        gene_count = np.array(genes_list)
        meanprophagelength = lengths.mean()
        meangenenumber = gene_count.mean()
        minprophagelength = lengths.min()
        mingenenumber = gene_count.min()
        maxprophagelength = lengths.max()
        maxgenenumber = gene_count.max()
        gc_lists = np.array(gc_list)
        meangc = gc_lists.mean()
        mingc = gc_lists.min()
        maxgc = gc_lists.max()
    else:
        meanprophagelength = 0
        meangenenumber = 0
        minprophagelength = 0
        mingenenumber = 0
        maxprophagelength = 0
        maxgenenumber = 0
        meangc = 0
        mingc = 0
        maxgc = 0

    # Add all information to results.tsv, for analysis, only want first two words in Organism
    with open('results.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow([f'{id}', f'{final}', f'{isolation_date}', f'{year}', f'{isolation_location}', f'{country}',
                             f'{isolation_source}', f'{prophage_num}', f'{meanprophagelength}', f'{minprophagelength}',
                             f'{maxprophagelength}',
                             f'{meangenenumber}', f'{mingenenumber}', f'{maxgenenumber}', f'{meangc}', f'{mingc}',
                             f'{maxgc}', f'{virulence_gene_prophages}', f'{virulence_genes}', f'{amr_gene_prophages}',
                             f'{amr_genes}'])
    os.chdir('..')
    subprocess.run(f"tar -zcf {id}.tar.gz {id}", shell=True)
    shutil.rmtree(id)
    print(f"{id} Done")

os.chdir('Output')
# Running mash on prophage nucleotide sequences to determine similarity to reference phages genomes
# os.chdir("fnasequences")
# subprocess.run("mash sketch -s 1000 -k 21 -o all_prophages *.fna", shell=True)
# subprocess.run("mv all_prophages.msh ..", shell=True)
# os.chdir("..")
# subprocess.run("mash dist -p 10 -v 0.001 -d 0.05 /data/DB/inphared5Jan/5Jan2023_genomes.fa.msh all_prophages.msh > distances.tab", shell=True)

if not type(args.start) is int:
    subprocess.run("tar -zcf fnasequences.tar.gz fnasequences", shell=True)
    subprocess.run("rm -r fnasequences", shell=True)
else:
    print(f'Start and end numbers given, if whole script ended, zip the fnasequences directory')

if args.skip_figures is False:
    # Create histogram showing all prophage lengths
    pylab.figure(figsize=(15, 9))  # Plot the count data
    pylab.hist(total_lengths, bins=20)
    pylab.title(
        "{0} Prophage Sequences\nLengths {1} to {2}".format(len(total_lengths), min(total_lengths), max(total_lengths))
    )
    pylab.xlabel("Sequence Length (bp)")
    pylab.ylabel("Count")
    pylab.rc('xtick', labelsize=12)
    pylab.rc('ytick', labelsize=12)
    pylab.savefig("prophagesequencelength")
    plt.close()

    # Create histogram showing all prophage gene counts
    pylab.figure(figsize=(15, 9))  # Plot the count data
    pylab.hist(total_genes, bins=20)
    pylab.title(
        "{0} Prophage Coding Sequence Counts\nCounts {1} to {2}".format(sum(total_genes), min(total_genes), max(total_genes))
    )
    pylab.xlabel("Number of Coding Sequences")
    pylab.ylabel("Count")
    pylab.rc('xtick', labelsize=12)
    pylab.rc('ytick', labelsize=12)
    pylab.savefig("prophagegenecounts")
    plt.close()

    # Create histogram showing prophages per genome counts
    data = pd.read_csv("results.tsv", sep="\t")
    prophage_number_counts = data['Number of Prophages']
    prophage_counts = prophage_number_counts.value_counts()
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(prophage_counts.index, prophage_counts.values, width=0.4)
    pylab.title("Number of Prophages per Genome")
    pylab.xlabel("Number of Prophages")
    pylab.ylabel("Count")
    pylab.rc('xtick', labelsize=12)
    pylab.rc('ytick', labelsize=12)
    pylab.xticks(range(min(prophage_counts.index), max(prophage_counts.index) + 1))
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("genomeprophagecounts")
    plt.close()

    # Create histogram showing virulence gene counts
    data = pd.read_csv("prophage_results.tsv", sep="\t")
    virulence_gene_counts = data['Virulence Gene Count']
    virulence_counts = virulence_gene_counts.value_counts()
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(virulence_counts.index, virulence_counts.values, width=0.4)
    pylab.title("Number of Virulence Genes per Prophage Genome")
    pylab.xlabel("Number of Virulence Genes")
    pylab.ylabel("Count")
    pylab.rc('xtick', labelsize=12)
    pylab.rc('ytick', labelsize=12)
    pylab.yscale('log')
    pylab.xticks(range(min(virulence_counts.index), max(virulence_counts.index) + 1))
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("prophagevirulencegenecounts")
    plt.close()

    # Create histogram showing amr gene counts
    data = pd.read_csv("prophage_results.tsv", sep="\t")
    amr_gene_counts = data['AMR Count']
    amr_counts = amr_gene_counts.value_counts()
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(amr_counts.index, amr_counts.values, width=0.4)
    pylab.title("Number of AMR genes per Prophage Genome")
    pylab.xlabel("Number of AMR Genes")
    pylab.ylabel("Count")
    pylab.rc('xtick', labelsize=12)
    pylab.rc('ytick', labelsize=12)
    pylab.yscale('log')
    pylab.xticks(range(min(amr_counts.index), max(amr_counts.index) + 1))
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("prophageamrgenecounts")
    plt.close()

    # Create graph for how many results there are for AMR genes against each drug class
    data = pd.read_csv("amr_genes.tsv", sep="\t")  # Read the tsv file into a DataFrame
    drug_classes = data['Drug Class']  # Extract the 'Drug Class' column
    class_counts = drug_classes.value_counts()  # Count the frequency of each drug class
    pylab.figure(figsize=(15, 9))  # Plot the count data
    pylab.tight_layout()
    bar_plot = pylab.bar(class_counts.index, class_counts.values)
    pylab.xlabel('Drug Class')
    pylab.ylabel('Count')
    pylab.title('Number of Results per Drug Class')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("prophagedrugclasscounts")
    plt.close()

    # Create graph for how many results there are for AMR genes for each cut-off (perfect, strict, loose)
    data = pd.read_csv("amr_genes.tsv", sep="\t")  # Read the tsv file into a DataFrame
    cut_off = data['Cut-Off']  # Extract the 'Drug Class' column
    cut_off_counts = cut_off.value_counts()  # Count the frequency of each drug class
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(cut_off_counts.index, cut_off_counts.values)
    pylab.xlabel('Cut-Off Level')
    pylab.ylabel('Count')
    pylab.title('Number of AMR Genes per Cut-Off Level')
    pylab.xticks(rotation=0)
    pylab.rc('xtick', labelsize=12)
    pylab.rc('ytick', labelsize=12)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("prophageamrcutoffcounts")
    plt.close()

    # Create histogram showing number of genomes for each species
    data = pd.read_csv("results.tsv", sep="\t")
    species = data['Species']
    species_counts = species.value_counts()
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(species_counts.index, species_counts.values)
    pylab.xlabel('Species Name')
    pylab.ylabel('Count')
    pylab.title('Number of Genomes per Enterococcus Species')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("speciescounts")
    plt.close()

    # Calculating and making bar charts for number of prophages per genome and prophage amr genes for each species
    results = pd.read_csv('results.tsv', sep='\t')
    grouped_data = results.groupby('Species')['Number of Prophages'].agg(['sum', 'mean']).reset_index()
    grouped_data.to_csv('SpeciesVProphages.csv', index=False)
    amr_grouped_data = results.groupby('Species')['Total number of AMR Genes in prophages'].agg(
        ['sum', 'mean']).reset_index()
    amr_grouped_data.to_csv('SpeciesVAMRs.csv', index=False)

    # Plot the bar chart for the sum of number of prophages for each species
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(grouped_data['Species'], grouped_data['sum'])
    pylab.xlabel('Species')
    pylab.ylabel('Sum of Number of Prophages')
    pylab.title('Sum of Number of Prophages per Species')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("prophagesperspeciescount")
    plt.close()

    # Plot the bar chart for the average number of prophages for each species
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(grouped_data['Species'], grouped_data['mean'])
    pylab.xlabel('Species')
    pylab.ylabel('Average Number of Prophages')
    pylab.title('Average Number of Prophages per Species')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.20)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        label = "{:.2f}".format(y_val)  # Format y_val to two decimal places
        pylab.text(x_val, y_val, label, ha='center', va='bottom')
    pylab.savefig("averageprophagesperspeciescount")
    plt.close()

    # Plot the bar chart for the sum of number of prophage AMR genes for each species
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(amr_grouped_data['Species'], amr_grouped_data['sum'])
    pylab.xlabel('Species')
    pylab.ylabel('Sum of Number of Prophage AMRs')
    pylab.title('Sum of Number of Prophage AMRs per Species')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("amrsperspeciescount")
    plt.close()

    # Plot the bar chart for the average number of prophage AMR genes for each species
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(amr_grouped_data['Species'], amr_grouped_data['mean'])
    pylab.xlabel('Species')
    pylab.ylabel('Average Number of Prophage AMRs per Genome')
    pylab.title('Average Number of Prophage AMRs per Genome by Species')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.20)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        label = "{:.2f}".format(y_val)  # Format y_val to two decimal places
        pylab.text(x_val, y_val, label, ha='center', va='bottom')
    pylab.savefig("averageamrperspeciescount")
    plt.close()

    # Create histogram showing number of genomes for each isolation Country
    data = pd.read_csv("results.tsv", sep="\t")
    country = data['Collection Country']
    country_counts = country.value_counts()
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(country_counts.index, country_counts.values)
    pylab.xlabel('Isolation Country')
    pylab.ylabel('Count')
    pylab.title('Number of Genomes per Isolation Country')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("countrycounts")
    plt.close()

    # Calculating and making bar charts for number of prophages per genome and prophage amr genes for each isolation country
    results = pd.read_csv('results.tsv', sep='\t')
    country_grouped_data = results.groupby('Collection Country')['Number of Prophages'].agg(['sum', 'mean']).reset_index()
    country_grouped_data.to_csv('CountryVProphages.csv', index=False)
    country_amr_grouped_data = results.groupby('Collection Country')['Total number of AMR Genes in prophages'].agg(
        ['sum', 'mean']).reset_index()
    country_amr_grouped_data.to_csv('CountryVAMRs.csv', index=False)

    # Plot the bar chart for the sum of number of prophages for each isolation country
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(country_grouped_data['Collection Country'], country_grouped_data['sum'])
    pylab.xlabel('Isolation Country')
    pylab.ylabel('Sum of Number of Prophages')
    pylab.title('Sum of Number of Prophages by Isolation Country')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("prophagesperisolationcountry")
    plt.close()

    # Plot the bar chart for the average number of prophages for each isolation country
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(country_grouped_data['Collection Country'], country_grouped_data['mean'])
    pylab.xlabel('Isolation Country')
    pylab.ylabel('Average Number of Prophages')
    pylab.title('Average Number of Prophages by Isolation Country')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.20)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        label = "{:.2f}".format(y_val)  # Format y_val to two decimal places
        pylab.text(x_val, y_val, label, ha='center', va='bottom')
    pylab.savefig("averageprophagesperisolationcountry")
    plt.close()

    # Plot the bar chart for the sum of number of prophage AMR genes for each isolation country
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(country_amr_grouped_data['Collection Country'], country_amr_grouped_data['sum'])
    pylab.xlabel('Isolation Country')
    pylab.ylabel('Sum of Number of Prophage AMRs')
    pylab.title('Sum of Number of Prophage AMRs by Isolation Country')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("amrsperisolationcountry")
    plt.close()

    # Plot the bar chart for the average number of prophage AMR genes for each isolation country
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(country_amr_grouped_data['Collection Country'], country_amr_grouped_data['mean'])
    pylab.xlabel('Isolation country')
    pylab.ylabel('Average Number of Prophage AMRs per Genome')
    pylab.title('Average Number of Prophage AMRs per Genome by Isolation Country')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.20)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        label = "{:.2f}".format(y_val)  # Format y_val to two decimal places
        pylab.text(x_val, y_val, label, ha='center', va='bottom')
    pylab.savefig("averageamrperisolationcountry")
    plt.close()

    # Create histogram showing number of genomes for each year of isolation
    data = pd.read_csv("results.tsv", sep="\t")
    year = data['Collection Year']
    year_counts = year.value_counts()
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(year_counts.index, year_counts.values)
    pylab.xlabel('Year of Isolation')
    pylab.ylabel('Count')
    pylab.title('Number of Genomes by Year of Isolation')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("yearcounts")
    plt.close()

    # Calculating and making bar charts for number of prophages per genome and prophage amr genes for each isolation year
    results = pd.read_csv('results.tsv', sep='\t')
    year_grouped_data = results.groupby('Collection Year')['Number of Prophages'].agg(['sum', 'mean']).reset_index()
    year_grouped_data.to_csv('YearVProphages.csv', index=False)
    year_amr_grouped_data = results.groupby('Collection Year')['Total number of AMR Genes in prophages'].agg(
        ['sum', 'mean']).reset_index()
    year_amr_grouped_data.to_csv('YearVAMRs.csv', index=False)

    # Plot the bar chart for the sum of number of prophages for each isolation year
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(year_grouped_data['Collection Year'], year_grouped_data['sum'])
    pylab.xlabel('Isolation Year')
    pylab.ylabel('Sum of Number of Prophages')
    pylab.title('Sum of Number of Prophages by Year of Isolation')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("prophagesperisolationyear")
    plt.close()

    # Plot the bar chart for the average number of prophages for each year of isolation
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(year_grouped_data['Collection Year'], year_grouped_data['mean'])
    pylab.xlabel('Isolation Year')
    pylab.ylabel('Average Number of Prophages')
    pylab.title('Average Number of Prophages by Isolation Year')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.20)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        label = "{:.2f}".format(y_val)  # Format y_val to two decimal places
        pylab.text(x_val, y_val, label, ha='center', va='bottom')
    pylab.savefig("averageprophagesperisolationyear")
    plt.close()

    # Plot the bar chart for the sum of number of prophage AMR genes for each year of isolation
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(year_amr_grouped_data['Collection Year'], year_amr_grouped_data['sum'])
    pylab.xlabel('Isolation Year')
    pylab.ylabel('Sum of Number of Prophage AMRs')
    pylab.title('Sum of Number of Prophage AMRs by Isolation Year')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.30)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
    pylab.savefig("amrsperisolationyear")
    plt.close()

    # Plot the bar chart for the average number of prophage AMR genes for each isolation year
    pylab.figure(figsize=(15, 9))  # Plot the count data
    bar_plot = pylab.bar(year_amr_grouped_data['Collection Year'], year_amr_grouped_data['mean'])
    pylab.xlabel('Isolation Year')
    pylab.ylabel('Average Number of Prophage AMRs per Genome')
    pylab.title('Average Number of Prophage AMRs per Genome by Isolation Year')
    pylab.xticks(rotation=20, ha='right')
    pylab.subplots_adjust(bottom=0.20)
    pylab.rc('xtick', labelsize=5)
    pylab.rc('ytick', labelsize=8)
    for bar in bar_plot:
        y_val = bar.get_height()
        x_val = bar.get_x() + bar.get_width() / 2
        label = "{:.2f}".format(y_val)  # Format y_val to two decimal places
        pylab.text(x_val, y_val, label, ha='center', va='bottom')
    pylab.savefig("averageleperisolationyear")
    plt.close()