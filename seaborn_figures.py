import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import matplotlib.ticker as ticker
from argparse import ArgumentParser
from matplotlib.ticker import LogLocator, FuncFormatter
import numpy as np


parser = argparse.ArgumentParser()

parser.add_argument("file", help="Input file(s)")
parser.add_argument("--xaxis", help="Column to use for x-axis data")
parser.add_argument("--yaxis", help="Column to use for y-axis data", required=False)
parser.add_argument("--Title", help="Graph Title")
parser.add_argument("--xlabel", help="Label for x-axis")
parser.add_argument("--ylabel", help="Label for y-axis")
parser.add_argument("--Box", help="Box Plot", action="store_true")
parser.add_argument("--Line", help="Line Graph", action="store_true")
parser.add_argument("--MultiLine", help="Line Graph with multiple groups to plot", action="store_true")
parser.add_argument("--Bar", help="Bar Chart", action="store_true")
parser.add_argument("--BarCount", help="Bar Chart for the Count of Specific numbers in column", action="store_true")
parser.add_argument("--Rel", help="Line Graph with regression", action="store_true")
parser.add_argument("--Violin", help="Violin Chart", action="store_true")
parser.add_argument("--Scatter", help="Scatter Graph", action="store_true")
parser.add_argument("--Pie", help="Pie Chart", action="store_true")
args = parser.parse_args()

name = args.file.split('.t')[0]

if args.Box:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)

    df_melted = pd.melt(df, id_vars="Phage", var_name="MOI", value_name="vals")

    plt.figure(figsize=[15, 9])
    sns.set_style("whitegrid")
    ax = sns.stripplot(data=df_melted, x="Phage", y="vals", hue="MOI",
                  size=10, linewidth=2, jitter=False, dodge=True, palette='bright')
    plt.ylabel(f'{args.ylabel}')
    plt.title(f'{args.Title}')
    plt.xticks(rotation=20, ha='center')
    plt.legend(title = "MOI", bbox_to_anchor = (1.02, 1),loc = "upper left",borderaxespad = 0)

    plt.tight_layout()
    #plt.show()
    plt.savefig(f"boxplot.png", format="png",bbox_inches="tight")
    plt.close()

if args.MultiLine:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)

    # Melt dataframe
    df_melted = pd.melt(df, id_vars="Time (minutes)", var_name="Antibiotic", value_name="vals")

    plt.figure(figsize=[15, 9])
    sns.set_style("whitegrid")
    ax = sns.pointplot(
        data=df_melted,
        x=f"{args.xaxis}",
        y="vals",
        hue="Antibiotic",
        palette="bright"
        )
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0
    )
    for collection in ax.collections:
        collection.set_sizes([5])
    plt.subplots_adjust(right=0.80)
    ax.set_ylim(bottom=0.01, top=10)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
    plt.yscale('log')
    plt.ylabel(f'{args.ylabel}')
    plt.xlabel(f'{args.xlabel}')
    plt.title(f'{args.Title}')
    plt.show()
    plt.close()


if args.Line:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)
    plt.figure(figsize=[15, 9])
    sns.set_style("whitegrid")
    ax = sns.scatterplot(
        data=df,
        x=f"{args.xaxis}", y=f"{args.yaxis}", palette="bright")
    ax.legend(fontsize='small')
    plt.subplots_adjust(right=0.70)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set_xlim(left=0, right=150)
    plt.yscale('log')
    plt.ylabel(f'{args.ylabel}')
    plt.xlabel(f'{args.xlabel}')
    plt.title(f'{args.Title}')
    plt.show()
    plt.close()

if args.Bar:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)
    sns.set_style("whitegrid")
    sns.barplot(data=df, x=f'{args.xaxis}', y=f'{args.yaxis}', color='steelblue')
    plt.yscale('log')
    plt.ylabel(f'{args.ylabel}')
    plt.xlabel(f'{args.xlabel}')
    plt.title(f'{args.Title}')
    plt.show()
    plt.close()

if args.BarCount:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)
    counts = df[f'{args.xaxis}'].value_counts().sort_index()
    sns.set_style("whitegrid")
    ax = sns.barplot(data=df, x=counts.index, y=counts.values, color='steelblue')
    plt.yscale('log')
    plt.ylabel(f'{args.ylabel}', fontsize=16, labelpad=15)
    plt.xlabel(f'{args.xlabel}', fontsize=16, labelpad=15)
    plt.title(f'{args.Title}', fontsize=16)
    # plt.xticks(rotation=20, ha='right', fontsize=12)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    for i, count in enumerate(counts.values):
        ax.text(i, count, str(count), ha='center', va='bottom', fontsize=14)
    plt.subplots_adjust(top=0.88, bottom=0.2)
    plt.show()
    plt.close()

if args.Rel:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)
    df.columns = (
        df.columns
        .str.replace(r'0\.', '0_', regex=True)  # Step 1: Replace "0." with "0_"
        .str.replace(r'\.\d+$', '', regex=True)  # Step 2: Remove any trailing ".1", ".2", etc.
        .str.replace(r'0_', '0.', regex=True)
    )

    # Melt the DataFrame to long format
    df_melted = pd.melt(df, id_vars=['Time (minutes)'], var_name='Group', value_name='Value')
    plt.figure(figsize=[10, 6])
    sns.set_style("whitegrid")
    sns.relplot(
        data=df_melted, kind="line", marker="o", linewidth=0,
        x="Time (minutes)", y="Value", hue="Group", palette="bright", errorbar="sd")
    plt.yscale('log')
    plt.ylabel('OD600 (log)')
    ax = plt.gca()
    ax.set_xlim(left=0, right=150)
    plt.yticks([0.1, 1, 10], ['0.1', '1', '10'])
    plt.savefig(f"{name}REG.svg", format="svg")
    #plt.show()
    plt.close()

if args.Violin:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)
    df_melted = df.melt(id_vars=['Genome Composition'],value_vars=['CpG Odds Ratio', 'ApT Odds Ratio', 'TpA Odds Ratio', 'GpC Odds Ratio',
                                    'ApA Odds Ratio', 'TpT Odds Ratio', 'ApC Odds Ratio',  'GpT Odds Ratio',
                                    'ApG Odds Ratio', 'CpT Odds Ratio', 'CpA Odds Ratio', 'TpG Odds Ratio',
                                    'CpC Odds Ratio', 'GpG Odds Ratio', 'GpA Odds Ratio', 'TpC Odds Ratio'],
                        var_name='Dinucleotide', value_name='Odds Ratio')
    df_melted['Dinucleotide'] = df_melted['Dinucleotide'].str.replace(' Odds Ratio', '')
    plt.figure(figsize=[15,9])
    sns.set_style("whitegrid")
    with sns.color_palette("bright"):
        if args.phage:
            sns.violinplot(data=df_melted, x='Dinucleotide', y='Odds Ratio', hue='Genome Composition', palette=["C4", "C8"])
        elif args.eukaryoticvirus:
            sns.violinplot(data=df_melted, x='Dinucleotide', y='Odds Ratio', hue='Genome Composition', palette=["C0", "C1", "C2", "C3"])
    plt.axhline(y=1.0, color='black', linewidth=1, linestyle='-', label="y=1.0")  # Solid line at y=1.0
    plt.axhline(y=0.75, color='red', linewidth=1, linestyle='--', label="y=0.75")  # Dashed line at y=0.75
    plt.axhline(y=1.25, color='blue', linewidth=1, linestyle='--', label="y=1.25")  # Dashed line at y=1.25
    plt.show()

if args.Scatter:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)

    #plt.figure(figsize=[10, 6])
    sns.set_theme()
    sns.lmplot(
        data=df,
        x=f"{args.xaxis}", y=f"{args.yaxis}",  scatter_kws={"s": 3, "alpha":0.3},
        line_kws={"ls": "--"}, palette="bright")
    plt.ylabel(f'{args.ylabel}')
    plt.xlabel(f'{args.xlabel}')
    plt.title(f'{args.Title}')
    plt.show()

if args.Pie:
    df = pd.read_csv(f"{args.file}", sep="\t", low_memory=False)

    value_counts = df[f'{args.xaxis}'].value_counts()
    # Get top 4
    top_4 = value_counts.nlargest(4)
    # Sum the rest
    other_sum = value_counts.iloc[4:].sum()
    pie_data = top_4.append(pd.Series({'Other': other_sum}))
    colours = sns.color_palette('bright', len(pie_data))
    fig, ax = plt.subplots(figsize=(14, 14))
    wedges, texts, autotexts = ax.pie(
        pie_data.values,
        labels=None,  # Don't show labels on wedges
        colors=colours,
        autopct='%.1f%%',
        pctdistance=0.9,  # Push percentages out
        startangle=0,
    )

    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(12)
    # Add legend with class names
    ax.legend(wedges, pie_data.index, title=f"{args.xaxis}", loc="center left", bbox_to_anchor=(1, 0.5), fontsize=14, title_fontsize=16)

    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()