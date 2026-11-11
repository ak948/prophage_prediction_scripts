import pandas as pd
import csv

with open('datatoadd.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['BiosampleID', 'Number of Prophages Carrying Endolysins', 'Prophage Endolysin Count',
                         'Number of Prophages Carrying AMGs', 'Prophage AMG Count',
                         'Number of Prophages Carrying Defence Systems', 'Prophage Defence System Count',
                         'Number of Prophages Carrying Anti-Defence Systems', 'Prophage Anti-Defence System Count'])

prophage_amgs = {}
prophage_endolysins = {}
prophage_antidefences = {}
prophage_defences = {}
prophages_with_amgs = {}
prophages_with_endolysins = {}
prophages_with_antidefences = {}
prophages_with_defences = {}

df = pd.read_csv('results.tsv', sep='\t')
biosampleid = df['BiosampleID']
for id in biosampleid:
    prophage_amgs[id] = 0
    prophage_endolysins[id] = 0
    prophage_antidefences[id] = 0
    prophage_defences[id] = 0
    prophages_with_amgs[id] = 0
    prophages_with_endolysins[id] = 0
    prophages_with_antidefences[id] = 0
    prophages_with_defences[id] = 0

with open('prophage_dataframe.tsv', 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        host_genome = row['BiosampleID']
        prophage_amgs[host_genome] += int(row['Prophage AMGs Count'])
        prophage_endolysins[host_genome] += int(row['Endolysin Count'])
        prophage_antidefences[host_genome] += int(row['Prophage Anti-Defense System Count'])
        prophage_defences[host_genome] += int(row['Prophage Defense System Count'])
        prophages_with_amgs[host_genome] += int(row['Prophage AMGs Present'])
        if int(row['Endolysin Count']) == 0:
            endolysin_present = 0
        else:
            endolysin_present = 1
        prophages_with_endolysins[host_genome] += endolysin_present
        prophages_with_antidefences[host_genome] += int(row['Prophage Anti-Defense System Present'])
        prophages_with_defences[host_genome] += int(row['Prophage Defense System Present'])

with open('datatoadd.tsv', 'a') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    for id in biosampleid:
        tsv_writer.writerow([f"{id}", prophages_with_endolysins[id], prophage_endolysins[id], prophages_with_amgs[id],
                             prophage_amgs[id], prophages_with_defences[id], prophage_defences[id],
                             prophages_with_antidefences[id], prophage_antidefences[id]])
