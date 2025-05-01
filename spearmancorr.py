import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import argparse
from argparse import ArgumentParser

parser = argparse.ArgumentParser()

parser.add_argument("--input", help="Input file(s)")
parser.add_argument("--Variable1", help="Name of Column in Input File containing First Variable")
parser.add_argument("--Variable2", help="Name of Column in Input File containing Second Variable")
parser.add_argument("--Output", help="Name of Output file prefix, i.e. Output_dunns.tsv and Output_tukey.tsv")
args = parser.parse_args()

df = pd.read_csv(f'{args.input}', sep='\t')

rho, p = spearmanr(df[args.Variable1], df[args.Variable2])
pear, _ = pearsonr(df[args.Variable1], df[args.Variable2])

print(f"rho = {rho}")
print(f"p = {p}")
print(f"pearson = {pear}")
if rho >0:
    print(f"Positive Correlation")
elif rho == 0:
    print(f"No Correlation")
elif rho < 0:
    print(f"Negative Correlation")
if p > 0.5:
    print(f"No Significant Difference")
elif p < 0.5:
    print(f"Significant Difference")