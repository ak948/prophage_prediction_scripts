import shutil
import argparse
from argparse import ArgumentParser
import pandas as pd
import os
import subprocess
import csv

parser = argparse.ArgumentParser()

parser.add_argument("input", help="path to input file containing list of BiosampleIDs")
parser.add_argument("accession", help="Name of column containing accession numbers in input file")
parser.add_argument("biosample", help="Name of column containing BiosampleIDs in input file")
parser.add_argument("cores", help="Number of Cores to use for VIBRANT and phageboost predictions", type=int)
parser.add_argument("--start", help="Number of BiosampleID to start on (For splitting into chunks)", type=int)
parser.add_argument("--end", help="Number of BiosampleID to end on (For splitting into chunks)", type=int)
parser.add_argument("email", help="Email to use for the Entrez.email function")
parser.add_argument("--phageboost", help="Location of phageboost executable, for if not in path", default="PhageBoost")
parser.add_argument("--vibrant", help="Location of VIBRANT executable, for if not in path", default="VIBRANT_run.py")
args = parser.parse_args()

if type(args.start) is int:
    df = pd.read_csv(args.input, sep='\t')
    assembly = df[args.accession].loc[args.start : args.end]
    biosampleid = df[args.biosample].loc[args.start : args.end]
else:
    df = pd.read_csv(args.input, sep='\t')
    assembly = df[args.accession]
    biosampleid = df[args.biosample]

if not os.path.isfile('genomeinfo.tsv'):
    with open('genomeinfo.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['Biosampleid', 'Organism', 'Submission Date', 'Last Update Date', 'Submitter', 'Organism Name', 'Collection Date', 'Collection Location', 'Isolation Source'])
if not os.path.isfile('repeats.txt'):
    with open('repeats.txt', 'a') as myfile:
        myfile.write('Biosampleid\tDate\tAccession\n')

from phageboost import *


for id, accession in zip(biosampleid, assembly):
    links = get_url(id, accession, args.email, download=True) #Get genome for the biosampleID
    assemblyinfo = get_assembly_info(id, args.email) #Get information on the actual sample the bacteria is from
    runphageboost = run_phageboost(id, args.cores, args.phageboost) #Run phageboost on the genome
    alsovibrant = also_vibrant(id, args.cores, args.vibrant) #Run vibrant on the prophage genomes predicted by phageboost
    os.chdir('..')
    subprocess.run("rm -r *.fna", shell=True) #Delete the genome that is no longer needed
    if os.path.isdir('phageboost_results'):
        os.chdir('..')
        subprocess.run(f"tar -zcf {id}.tar.gz {id}", shell=True)
        shutil.rmtree(id)
    elif os.path.isdir('prophage_prediction'):
        os.chdir('prophage_prediction')


print("Script finished, change to environment containing genomad before changing to next script")

