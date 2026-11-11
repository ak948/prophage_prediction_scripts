import pandas as pd
import csv

"""
Read vOTUs_with_multispecie_hostrange.tsv to get list of vOTUs with host range of multiple species
For each prophage in prophage_dataframe.tsv, if 'Enterococcus Prophage Species' matches one of the vOTUs, add BiosampleID to dict after checking if already there
Chuck list of BiosampleIDs in new spreadsheet? 
"""

dt = pd.read_csv(f'vOTUs_with_multispecie_hostrange.tsv', sep='\t')
multispeciehost_votus = set(dt['vOTU'].astype(str).str.strip())
prophage_with_multispecies = {}

with open('../prophage_dataframe.tsv', 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        votu = row['Enterococcus Prophage Species'].strip()
        if votu in multispeciehost_votus:
            genome = row['BiosampleID']
            if genome in prophage_with_multispecies:
                if prophage_with_multispecies[genome] == '0':
                    prophage_with_multispecies[genome] = '1'
            else:
                prophage_with_multispecies[genome] = '1'
        else:
            genome = row['BiosampleID']
            if genome not in prophage_with_multispecies:
                prophage_with_multispecies[genome] = '0'

with open('genomes_with_prophages_infecting_multispecies.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['BiosampleID', 'Prophage Infecting Multiple Species'])  # Header
    for key, value in prophage_with_multispecies.items():
        writer.writerow([key, value])
