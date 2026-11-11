import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = argparse.ArgumentParser()

parser.add_argument("--input", help="Input file(s)")
parser.add_argument("--Variable1", help="Name of Column in Input File containing First Variable")
parser.add_argument("--Variable2", help="Name of Column in Input File containing Second Variable")
parser.add_argument("--Output", help="Name of Output file prefix, i.e. Output_dunns.tsv and Output_tukey.tsv")
args = parser.parse_args()

df = pd.read_csv(f'{args.input}', sep='\t', low_memory=False)

rho, p = spearmanr(df[args.Variable1], df[args.Variable2])
spearman_corr = spearmanr(df[args.Variable1], df[args.Variable2])
sns.heatmap(spearman_corr, annot=True, comap='viridis')
plot.title(f"Spearman Correlation Heatmap - {args.Variable1} vs. {args.Variable2}")
plt.show()

pear, _ = pearsonr(df[args.Variable1], df[args.Variable2])
pearson_corr = pearsonr(df[args.Variable1], df[args.Variable2])
sns.heatmap(pearson_corr, annot=True, comap='coolwarm')
plot.title(f"Pearson Correlation Heatmap - {args.Variable1} vs. {args.Variable2}")
plt.show()

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