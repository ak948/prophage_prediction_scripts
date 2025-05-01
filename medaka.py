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
parser.add_argument("--model", help="Nanopore model used", default="r1041_e82_400bps_fast_g615")
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t')
barcodes = df[args.barcodes]
samplenames = df[args.samples]

for barcode, sample in zip(barcodes, samplenames):
    subprocess.run(f'cd {sample}', shell=True)
    subprocess.run(f'medaka_consensus -i input_reads.fastq -d flye_{sample}/assembly.fasta -o polish1 -m {args.model}',
                   shell=True)
    subprocess.run(f'medaka_consensus -i input_reads.fastq -d polish1/consensus.fasta -o polish2 -m {args.model}',
                   shell=True)
    subprocess.run(f'medaka_consensus -i input_reads.fastq -d polish2/consensus.fasta -o polish3 -m {args.model}',
                   shell=True)
    subprocess.run(f'medaka_consensus -i input_reads.fastq -d polish3/consensus.fasta -o polish4 -m {args.model}',
                   shell=True)
    subprocess.run(f'cp polish4/consensus.fasta {sample}.fasta && gzip input_reads.fastq', shell=True)
    subprocess.run(f'rm -r flye_{sample}/00* flye_{sample}/10* flye_{sample}/20* flye_{sample}/30* flye_{sample}/40* '
                   f'flye_{sample}/assembly.fasta.* flye_{sample}/assembly_* flye_{sample}/params.json polish*/*.bam* '
                   f'polish*/consensus.fasta.* polish*/consensus_probs.hdf', shell=True)
    subprocess.run(f'cd ..', shell=True)