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


with open('gtdbresults.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['ID', 'Classification'])

with open('gtdbbac.tsv', 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        # Extract the desired data from different columns
        id = row['user_genome']
        species = row['classification']
        # Splitting the organism string to only keep the species
        split_species = species.split(";s__")

        # Extracting the species
        classification = split_species[1]

        if len(classification) > 0:
            with open('gtdbresults.tsv', 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow([f'{id}', f'{classification}'])
        else:
            split_species = species.split(";g__")

            # Extracting the species
            classification = split_species[1]
            with open('gtdbresults.tsv', 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow([f'{id}', f'{classification}'])