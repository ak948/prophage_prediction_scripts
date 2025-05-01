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
import sys

with open('../results.tsv', 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        name = row['Biosampleid']
        max_prophage = row['Number of Prophages']
        prophage_num = 1
        while int(prophage_num) <= int(max_prophage):
            subprocess.run(f"sed -i 's/gnl|Prokka|{name}_p{prophage_num}_1/{name}_p{prophage_num}/g' {name}_p{prophage_num}.fna", shell=True)
            prophage_num += 1
        print(f'{name} done')
