"""
Open gbk file, for each gene in it, get the gene name and add to count for that specific gene name, save as current_gene
and move previous_gene from current_gene to previous_gene (Do in reverse so is kept). If current_gene is endolysin,
save info to separate spreadsheet along with previous and next gene, and the sequence
"""

from Bio import SeqIO
from collections import defaultdict
import os
import csv


gene_counts = defaultdict(int)

with open('prophage_endolysins.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['Endolysin Name', 'Prophage Name', 'BiosampleID', 'Gene Before', 'Gene After',
                         'Endolysin Sequence'])
with open('phrogs_genes.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['PHROGs Gene', 'Count'])

os.chdir('FinalProphagesGbk')
for file in os.listdir():
    if file.endswith('.gbk'):
        prophage_name = file.split('.')[0]
        biosampleid = prophage_name.split('_')[0]
        endolysin_count = 0
        gene_count = 0
        previous_gene = ""
        current_gene = ""

        for seq_record in SeqIO.parse(f"{prophage_name}.gbk", "genbank"):
            for feature in seq_record.features:
                next_gene = ""
                if "product" in feature.qualifiers:
                    previous_gene = current_gene
                    current_gene = feature.qualifiers["product"][0]
                    gene_counts[current_gene] += 1
                    if current_gene == "endolysin":
                        endolysin_count += 1
                        try:
                            next_gene = seq_record.features[gene_count + 1].qualifiers["product"][0]
                        except:
                            continue
                        endolysin_sequence = seq_record.features[gene_count].qualifiers["translation"][0]
                        with open('../prophage_endolysins.tsv', 'a') as outfile:
                            tsv_writer = csv.writer(outfile, delimiter='\t')
                            tsv_writer.writerow(
                                [f'{prophage_name}_{endolysin_count}', f'{prophage_name}', f'{biosampleid}',
                                 f'{previous_gene}', f'{next_gene}', f'{endolysin_sequence}'])
                gene_count += 1

with open('../phrogs_genes.tsv', 'a') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    for gene in gene_counts:
        tsv_writer.writerow([f"{gene}", gene_counts[gene]])
