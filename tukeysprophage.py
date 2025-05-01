import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp
from sklearn.datasets import load_iris
import argparse
from argparse import ArgumentParser

parser = argparse.ArgumentParser()

parser.add_argument("--input", help="Input file(s)")
parser.add_argument("--Group", help="Name of Column in Input File containing Groups")
parser.add_argument("--Variable", help="Name of Column in Input File containing Variable being looked at for significance")
parser.add_argument("--Output", help="Name of Output file prefix, i.e. Output_dunns.tsv and Output_tukey.tsv")
args = parser.parse_args()

df = pd.read_csv(f'{args.input}', sep='\t')

group_counts = df[f'{args.Group}'].value_counts()
valid_groups = group_counts[group_counts >= 10].index

df_filtered = df[df[f'{args.Group}'].isin(valid_groups)]
targets_df = df_filtered[[f'{args.Group}',f'{args.Variable}']].copy()

grouped_data = df_filtered.groupby(f'{args.Group}')[f'{args.Variable}'].apply(list)

# Perform the one-way ANOVA
f_statistic, p_value = f_oneway(*grouped_data)

# Print the results
print('F-statistic:', f_statistic)
print('p-value:', p_value)

tukey = pairwise_tukeyhsd(endog=df_filtered[f'{args.Variable}'],
                          groups=df_filtered[f'{args.Group}'],
                          alpha=0.05)

print(tukey)

dg = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
dg.to_csv(f'{args.Output}_tukey.tsv', sep='\t', index=False)

dunn_and_done = sp.posthoc_dunn(targets_df, val_col = f'{args.Variable}', group_col = f'{args.Group}', p_adjust='holm')

dunn_and_done.to_csv(f'{args.Output}_dunns.tsv', sep='\t', index=False)