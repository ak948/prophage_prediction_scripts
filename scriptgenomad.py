import shutil
import argparse
from argparse import ArgumentParser
import pandas as pd
import os
import subprocess
import csv

parser = argparse.ArgumentParser()

parser.add_argument("--genomaddb", help="Location of genomad database, for if not in path")
parser.add_argument("prokkadb", help="Location of prokka database, automatically set to use a HMM database")
parser.add_argument("cores", help="Number of Cores to use for VIBRANT and phageboost predictions", type=int)
parser.add_argument("--start", help="Number of BiosampleID to start on (For splitting into chunks)", type=int)
parser.add_argument("--end", help="Number of BiosampleID to end on (For splitting into chunks)", type=int)
parser.add_argument("--genomad", help="Location of genomad executable, for if not in path", default="genomad")
parser.add_argument("--prokka", help="Location of prokka executable, for if not in path", default="prokka")
args = parser.parse_args()

if type(args.start) is int:
    df = pd.read_csv('genomeinfo.tsv', sep='\t')
    biosampleid = df['Biosampleid'].loc[args.start : args.end]
else:
    df = pd.read_csv('genomeinfo.tsv', sep='\t')
    biosampleid = df['Biosampleid']

from genomad import *

for id in biosampleid:
    subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)
    os.remove(f'{id}.tar.gz')
    rungenomad = run_genomad(id, args.genomaddb, args.cores, args.genomad)
    keepgenomes = predicted_prophages(id)
    runprokka = run_prokka(id, args.prokkadb, args.cores, args.prokka)
    subprocess.run("rm -r final*", shell=True)
    os.chdir('..')
    os.chdir('..')
    subprocess.run(f"tar -zcf {id}.tar.gz {id}", shell=True)
    shutil.rmtree(id)
