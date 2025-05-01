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


# os.chdir('Completed')
# Get list of biosampleIDs
df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid'].loc[0:1000]
# Write headings for file results for each genome will be recorded in
with open('results.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(
        ['Biosampleid', 'Species', 'Collection Date', 'Collection Location', 'Isolation Source', 'Number of Prophages',
         'Average Prophage Length', 'Minimum Prophage Length', 'Maximum Prophage Length',
         'Average number of prophage coding sequences', 'Minimum number of prophage coding sequences',
         'Maximum number of prophage coding sequences', 'Number of prophages containing AMR Genes',
         'Total number of AMR Genes in prophages'])
# Write headings for file the rgi-card results will be recorded in for each genome
with open('amr_genes.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['AMR Gene ID', 'ORF ID', 'Source Sequence', 'Cut-Off', 'Drug Class', 'Resistance Mechanism',
                         'AMR Gene Family', '% Length of Reference Genome'])

total_lengths = []
total_genes = []
total_amr_genes = []

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
                isolation_source = row['Isolation Source']
                isolation_location = row['Collection Location']

    os.chdir(f'{id}/phageboost_results')

    prophage_num = 0
    amr_gene_prophages = 0
    amr_genes = 0
    length_list = []
    genes_list = []
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
            # Run genome through RGI-CARD
            subprocess.run(f"rgi main -i {id}_p{file_num}.fna -o {id}_p{file_num}_rgi --clean", shell=True)
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
                    drug_class = row['Drug Class']
                    resistance = row['Resistance Mechanism']
                    amr_gene_family = row['AMR Gene Family']
                    percentage_length = row['Percentage Length of Reference Sequence']
                    os.chdir('../../..')
                    with open('amr_genes.tsv', 'a') as out_file:
                        tsv_writer = csv.writer(out_file, delimiter='\t')
                        tsv_writer.writerow([f'{amr_gene_id}', f'{orf_id}', f'{source_sequence}', f'{cut_off}',
                                             f'{drug_class}', f'{resistance}', f'{amr_gene_family}',
                                             f'{percentage_length}'])
                    row_num += 1
                    amr_genes += 1
                    prophage_amr_genes += 1
                    os.chdir(f'{id}/phageboost_results/{id}_p{file_num}')

            os.chdir('..')
            file_num += 1
            if prophage_amr_genes > 0:
                amr_gene_prophages += 1
        total_amr_genes.append(int(amr_genes))

    os.chdir('../..')
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
    else:
        meanprophagelength = 0
        meangenenumber = 0
        minprophagelength = 0
        mingenenumber = 0
        maxprophagelength = 0
        maxgenenumber = 0

    # Add all information to results.tsv, for analysis, only want first two words in Organism
    with open('results.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow([f'{id}', f'{final}', f'{isolation_date}', f'{isolation_location}', f'{isolation_source}',
                             f'{prophage_num}', f'{meanprophagelength}', f'{minprophagelength}', f'{maxprophagelength}',
                             f'{meangenenumber}', f'{mingenenumber}', f'{maxgenenumber}',
                             f'{amr_gene_prophages}', f'{amr_genes}'])
    subprocess.run(f"tar -zcf {id}.tar.gz {id}", shell=True)
    shutil.rmtree(id)
    print(f"{id} Done")

# Create histogram showing all prophage lengths
pylab.figure(figsize=(15, 9))  # Plot the count data
pylab.hist(total_lengths, bins=20)
pylab.title(
    "%i prophage sequences\nLengths %i to %i" % (len(total_lengths), min(total_lengths), max(total_lengths))
)
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.rc('xtick', labelsize=12)
pylab.rc('ytick', labelsize=12)
pylab.savefig("prophagesequencelength")

# Create histogram showing all prophage gene counts
pylab.figure(figsize=(15, 9))  # Plot the count data
pylab.hist(total_genes, bins=20)
pylab.title(
    "%i prophage coding sequence counts\nCounts %i to %i" % (sum(total_genes), min(total_genes), max(total_genes))
)
pylab.xlabel("Number of Coding Sequences")
pylab.ylabel("Count")
pylab.rc('xtick', labelsize=12)
pylab.rc('ytick', labelsize=12)
pylab.savefig("prophagegenecounts")

# Create histogram showing amr gene counts
data = pd.read_csv("results.tsv", sep="\t")
amr_gene_counts = data['Total number of AMR Genes in prophages']
amr_counts = amr_gene_counts.value_counts()
pylab.figure(figsize=(15, 9))  # Plot the count data
bar_plot = pylab.bar(amr_counts.index, amr_counts.values, width=0.4)
pylab.title("Number of amr genes per genome")
pylab.xlabel("Number of AMR Genes")
pylab.ylabel("Count")
pylab.rc('xtick', labelsize=12)
pylab.rc('ytick', labelsize=12)
pylab.xticks(range(min(amr_counts.index), max(amr_counts.index)+1))
for bar in bar_plot:
    y_val = bar.get_height()
    x_val = bar.get_x() + bar.get_width() / 2
    pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
pylab.savefig("prophageamrgenecounts")

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

# Create graph for how many results there are for AMR genes for each cut-off (perfect, strict, loose)
data = pd.read_csv("amr_genes.tsv", sep="\t")  # Read the tsv file into a DataFrame
drug_classes = data['Cut-Off']  # Extract the 'Drug Class' column
class_counts = drug_classes.value_counts()  # Count the frequency of each drug class
pylab.figure(figsize=(15, 9))  # Plot the count data
bar_plot = pylab.bar(class_counts.index, class_counts.values)
pylab.xlabel('Cut-Off Class')
pylab.ylabel('Count')
pylab.title('Number of AMR Genes for each cut-off level')
pylab.xticks(rotation=0)
pylab.rc('xtick', labelsize=12)
pylab.rc('ytick', labelsize=12)
for bar in bar_plot:
    y_val = bar.get_height()
    x_val = bar.get_x() + bar.get_width() / 2
    pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
pylab.savefig("prophageamrcutoffcounts")

# Create histogram showing number of genomes for each species
data = pd.read_csv("results.tsv", sep="\t")
species = data['Species']
species_counts = species.value_counts()
pylab.figure(figsize=(15, 9))  # Plot the count data
bar_plot = pylab.bar(species_counts.index, species_counts.values)
pylab.xlabel('Species Name')
pylab.ylabel('Count')
pylab.title('Number of genomes for each Enterococcus species')
pylab.xticks(rotation=20, ha='right')
pylab.subplots_adjust(bottom=0.30)
pylab.rc('xtick', labelsize=5)
pylab.rc('ytick', labelsize=8)
for bar in bar_plot:
    y_val = bar.get_height()
    x_val = bar.get_x() + bar.get_width() / 2
    pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
pylab.savefig("speciescounts")

# Calculating and making bar charts for number of prophages per genome for each species
results = pd.read_csv('results.tsv', sep='\t')
grouped_data = results.groupby('Species')['Number of Prophages'].agg(['sum', 'mean']).reset_index()

# Plot the bar chart for the sum of number of prophages
pylab.figure(figsize=(15, 9))  # Plot the count data
bar_plot = pylab.bar(grouped_data['Species'], grouped_data['sum'])
pylab.xlabel('Species')
pylab.ylabel('Sum of Number of Prophages')
pylab.title('Sum of Number of Prophages by Species')
pylab.xticks(rotation=20, ha='right')
pylab.subplots_adjust(bottom=0.30)
pylab.rc('xtick', labelsize=5)
pylab.rc('ytick', labelsize=8)
for bar in bar_plot:
    y_val = bar.get_height()
    x_val = bar.get_x() + bar.get_width() / 2
    pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom')
pylab.savefig("prophagesperspeciescount")

# Plot the bar chart for the average number of prophages
pylab.figure(figsize=(15, 9))  # Plot the count data
bar_plot = pylab.bar(grouped_data['Species'], grouped_data['mean'])
pylab.xlabel('Species')
pylab.ylabel('Average Number of Prophages')
pylab.title('Average Number of Prophages by Species')
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