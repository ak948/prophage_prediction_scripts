import matplotlib.pyplot as plt
import pylab
import numpy as np
import pandas as pd
import os
import csv

# Create histogram showing virulence gene counts
data = pd.read_csv("prophage_dataframe.tsv", sep="\t")
virulence_gene_counts = data['Prophage Virulence Gene Count']
virulence_counts = virulence_gene_counts.value_counts()
pylab.figure(figsize=(15, 9))  # Plot the count data
bar_plot = pylab.bar(virulence_counts.index, virulence_counts.values, width=0.4)
pylab.yscale('log')
pylab.title("Number of Virulence Genes per Prophage Genome", fontsize=12)
pylab.xlabel("Number of Virulence Genes", fontsize=12, labelpad=20)
pylab.ylabel("Count", fontsize=12)
pylab.rc('xtick', labelsize=16)
pylab.rc('ytick', labelsize=16)
pylab.xticks(range(min(virulence_counts.index), max(virulence_counts.index) + 1))
for bar in bar_plot:
    y_val = bar.get_height()
    x_val = bar.get_x() + bar.get_width() / 2
    pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom', fontsize=12)
pylab.savefig("prophagevirulencegenecounts")
plt.close()

# Create histogram showing amr gene counts
data = pd.read_csv("prophage_dataframe.tsv", sep="\t")
amr_gene_counts = data['Prophage ARG Count']
amr_counts = amr_gene_counts.value_counts()
pylab.figure(figsize=(15, 9))  # Plot the count data
bar_plot = pylab.bar(amr_counts.index, amr_counts.values, width=0.4)
pylab.yscale('log')
pylab.title("Number of AMR genes per Prophage Genome", fontsize=12)
pylab.xlabel("Number of AMR Genes", fontsize=13, labelpad=20)
pylab.ylabel("Count", fontsize=14)
pylab.rc('xtick', labelsize=14)
pylab.rc('ytick', labelsize=14)
pylab.xticks(range(min(amr_counts.index), max(amr_counts.index) + 1))
for bar in bar_plot:
    y_val = bar.get_height()
    x_val = bar.get_x() + bar.get_width() / 2
    pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom', fontsize=12)
pylab.savefig("prophageamrgenecounts")
plt.close()

data = pd.read_csv("genus.tsv", sep="\t")  # Read the tsv file into a DataFrame
taxonomy = data['Prophage Taxonomy']  # Extract the 'Drug Class' column
taxonomy_counts = data['Count']  # Count the frequency of each drug class
pylab.figure(figsize=(15, 9))  # Plot the count data
pylab.tight_layout()
bar_plot = pylab.bar(taxonomy, taxonomy_counts, width=0.4)
pylab.yscale('log')
pylab.xlabel('Prophage Taxonomy', fontsize=14, labelpad=20)
pylab.ylabel('Count', fontsize=15)
pylab.title('Prophage Taxonomy', fontsize=12)
pylab.xticks(rotation=0, ha='center')
pylab.rc('xtick', labelsize=12)
pylab.rc('ytick', labelsize=15)
for bar in bar_plot:
    y_val = bar.get_height()
    x_val = bar.get_x() + bar.get_width() / 2
    pylab.text(x_val, y_val, str(y_val), ha='center', va='bottom', fontsize=12)
pylab.savefig("prophagetaxonomycounts")
plt.close()
