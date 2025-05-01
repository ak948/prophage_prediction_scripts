
import csv
import os
import pandas as pd
import argparse
from argparse import ArgumentParser

parser = argparse.ArgumentParser()

parser.add_argument("input", help="Similarity table file")
parser.add_argument("program", help="Was similarity done with Maniac or taxmyphage", choices=('maniac', 'taxmyphage'))
args = parser.parse_args()

with open('taxonomy_similarity_clusters.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(["Prophage Name", "Phage Genus", "Phage Species"])

genus_dict = {}
species_dict = {}

if args.program == 'taxmyphage':
    with open(f'{args.input}', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        row_num = 1
        phage_genus = 1
        phage_species = 1
        for row in reader:
            row_num += 1
            sampleA = row['A']
            sampleB = row['B']
            similarity = float(row['sim'])
            if sampleB not in genus_dict:
                genus_dict[sampleB] = {}
            if sampleA not in genus_dict:
                genus_dict[sampleA] = {}
                if float(similarity) >= 70:
                    if genus_dict[sampleB]:
                        genus = genus_dict[sampleB]
                        genus_dict[sampleA] = genus
                    else:
                        genus_dict[sampleA] = phage_genus
                        phage_genus += 1
                        genus = genus_dict[sampleA]
                        genus_dict[sampleB] = genus
                else:
                    genus_dict[sampleA] = phage_genus
                    phage_genus += 1

                with open(f'{args.input}', 'r') as file3:
                    reader3 = csv.DictReader(file3, delimiter='\t')
                    row_number = 1
                    for row3 in reader3:
                        row_number +=1
                        if row_number >= row_num:
                            if row3['A'] == sampleA:
                                if float(row3['sim']) >= 70:
                                    sampleC = row3['B']
                                    if sampleC not in genus_dict:
                                        genus_dict[sampleC] = {}
                                    if not genus_dict[sampleC]:
                                        genus = genus_dict[sampleA]
                                        genus_dict[sampleC] = genus
                            elif row3['B'] == sampleA:
                                if float(row3['sim']) >= 70:
                                    sampleC = row3['A']
                                    if sampleC not in genus_dict:
                                        genus_dict[sampleC] = {}
                                    if not genus_dict[sampleC]:
                                        genus = genus_dict[sampleA]
                                        genus_dict[sampleC] = genus

            if sampleB not in species_dict:
                species_dict[sampleB] = {}
            if sampleA not in species_dict:
                species_dict[sampleA] = {}
                if float(similarity) >= 95:
                    if species_dict[sampleB]:
                        species = species_dict[sampleB]
                        species_dict[sampleA] = species
                    else:
                        species_dict[sampleA] = phage_species
                        phage_species += 1
                        species = species_dict[sampleA]
                        species_dict[sampleB] = species
                else:
                    species_dict[sampleA] = phage_species
                    phage_species += 1

                with open(f'{args.input}', 'r') as file3:
                    reader3 = csv.DictReader(file3, delimiter='\t')
                    row_number = 1
                    for row3 in reader3:
                        row_number += 1
                        if row_number >= row_num:
                            if row3['A'] == sampleA:
                                if float(row3['sim']) >= 95:
                                    sampleC = row3['B']
                                    if sampleC not in species_dict:
                                        species_dict[sampleC] = {}
                                    if not species_dict[sampleC]:
                                        species = species_dict[sampleA]
                                        species_dict[sampleC] = species
                            elif row3['B'] == sampleA:
                                if float(row3['sim']) >= 95:
                                    sampleC = row3['A']
                                    if sampleC not in species_dict:
                                        species_dict[sampleC] = {}
                                    if not species_dict[sampleC]:
                                        species = species_dict[sampleA]
                                        species_dict[sampleC] = species
elif args.program == 'maniac':
    with open(f'{args.input}', 'r') as file:
        reader = csv.DictReader(file, delimiter=',')
        row_num = 1
        phage_genus = 1
        phage_species = 1
        for row in reader:
            row_num =+1
            sampleA = row['Seq1']
            sampleB = row['Seq2']
            similarity = float(row['ANI'])
            if sampleB not in genus_dict:
                genus_dict[sampleB] = {}
            if sampleA not in genus_dict:
                genus_dict[sampleA] = {}
                if float(similarity) >= 0.70:
                    if genus_dict[sampleB]:
                        genus = genus_dict[sampleB]
                        genus_dict[sampleA] = genus
                    else:
                        genus_dict[sampleA] = phage_genus
                        phage_genus += 1
                        genus = genus_dict[sampleA]
                        genus_dict[sampleB] = genus
                else:
                    genus_dict[sampleA] = phage_genus
                    phage_genus += 1

                with open(f'{args.input}', 'r') as file3:
                    reader3 = csv.DictReader(file3, delimiter=',')
                    row_number = 1
                    for row3 in reader3:
                        row_number +=1
                        if row_number >= row_num:
                            if row3['Seq1'] == sampleA:
                                if float(row3['ANI']) >= 0.70:
                                    sampleC = row3['Seq2']
                                    if sampleC not in genus_dict:
                                        genus_dict[sampleC] = {}
                                    if not genus_dict[sampleC]:
                                        genus = genus_dict[sampleA]
                                        genus_dict[sampleC] = genus
                            elif row3['Seq2'] == sampleA:
                                if float(row3['ANI']) >= 0.70:
                                    sampleC = row3['Seq1']
                                    if sampleC not in genus_dict:
                                        genus_dict[sampleC] = {}
                                    if not genus_dict[sampleC]:
                                        genus = genus_dict[sampleA]
                                        genus_dict[sampleC] = genus

            if sampleB not in species_dict:
                species_dict[sampleB] = {}
            if sampleA not in species_dict:
                species_dict[sampleA] = {}
                if float(similarity) >= 0.95:
                    if species_dict[sampleB]:
                        species = species_dict[sampleB]
                        species_dict[sampleA] = species
                    else:
                        species_dict[sampleA] = phage_species
                        phage_species += 1
                        species = species_dict[sampleA]
                        species_dict[sampleB] = species
                else:
                    species_dict[sampleA] = phage_species
                    phage_species += 1

                with open(f'{args.input}', 'r') as file3:
                    reader3 = csv.DictReader(file3, delimiter=',')
                    row_number = 1
                    for row3 in reader3:
                        row_number +=1
                        if row_number >= row_num:
                            if row3['Seq1'] == sampleA:
                                if float(row3['ANI']) >= 0.95:
                                    sampleC = row3['Seq2']
                                    if sampleC not in species_dict:
                                        species_dict[sampleC] = {}
                                    if not species_dict[sampleC]:
                                        species = species_dict[sampleA]
                                        species_dict[sampleC] = species
                            elif row3['Seq2'] == sampleA:
                                if float(row3['ANI']) >= 0.95:
                                    sampleC = row3['Seq1']
                                    if sampleC not in species_dict:
                                        species_dict[sampleC] = {}
                                    if not species_dict[sampleC]:
                                        species = species_dict[sampleA]
                                        species_dict[sampleC] = species


with open('taxonomy_similarity_clusters.tsv', 'a') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    for sample in genus_dict:
        if not genus_dict[sample]:
            genus_dict[sample] = phage_genus
            phage_genus += 1
        if not species_dict[sample]:
            species_dict[sample] = phage_species
            phage_species += 1
        tsv_writer.writerow([sample, genus_dict[sample], species_dict[sample]])
