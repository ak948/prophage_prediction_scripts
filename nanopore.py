import argparse
from argparse import ArgumentParser
import pandas as pd
import os
import subprocess
import csv

parser = argparse.ArgumentParser()

parser.add_argument("input", help="path to input file containing list of barcodes")
parser.add_argument("--barcodes", help="Name of column containing barcode names", default="Barcodes")
parser.add_argument("--samples", help="Name of column containing sample names", default="Sample Names")
parser.add_argument("directory", help="Name of directory to put files in")
parser.add_argument("threads", help="Number of threads to use")
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t')
barcodes = df[args.barcodes]
samplenames = df[args.samples]

os.mkdir(args.directory)
subprocess.run(f'mv {input} *_fastq.gz {args.directory} && cd {args.directory}', shell=True)

for barcode, sample in zip(barcodes, samplenames):
    os.mkdir(sample)
    subprocess.run(f'mv *_{barcode}.fastq.gz {sample}', shell=True)
    subprocess.run(f'cd {sample}', shell=True)
    subprocess.run(f'filtlong --min_length 1000 --keep_percent 90 --target_bases 10000000 *_{barcode}.fastq.gz | gzip >'
                   f' output.fastq.gz && gunzip output.fastq.gz', shell=True)
    subprocess.run(f'flye --nano-raw output.fastq --threads {args.threads} --out-dir flye_{sample} && gzip output.fastq'
                   f' && gunzip input_reads.fastq.gz', shell=True)
    subprocess.run(f'cd ..', shell=True)
