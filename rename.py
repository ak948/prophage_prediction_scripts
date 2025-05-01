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

with open('Output/results.tsv', 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        name = row['Biosampleid']
        max_prophage = row['Number of Prophages']
        prophage_num = 1
        subprocess.run(f'tar -zxf {name}.tar.gz {name}', shell=True)
        subprocess.run(f'rm {name}.tar.gz', shell=True)
        os.chdir(f'{name}/phageboost_results')
        while prophage_num <= max_prophage:
            os.chdir(f'{name}_p{prophage_num}')
            subprocess.run(f"sed -i 's/gnl|Prokka|{name}_p{prophage_num}_1/{name}_p{prophage_num}/g' {name}_p{prophage_num}.fna", shell=True)
            subprocess.run(
                f"sed -i 's/gnl|Prokka|{name}_p{prophage_num}_1/{name}_p{prophage_num}/g' {name}_p{prophage_num}.faa",
                shell=True)
            subprocess.run(
                f"sed -i 's/gnl|Prokka|{name}_p{prophage_num}_1/{name}_p{prophage_num}/g' {name}_p{prophage_num}.gbk",
                shell=True)
            subprocess.run(
                f"sed -i 's/gnl|Prokka|{name}_p{prophage_num}_1/{name}_p{prophage_num}/g' {name}_p{prophage_num}.gff",
                shell=True)
            os.chdir(f'..')
            prophage_num += 1
        os.chdir('../..')
        subprocess.run(f'tar -zcf {name}.tar.gz {name}', shell=True)
        subprocess.run(f'rm -r {name}', shell=True)
        print(f'{name} done')
