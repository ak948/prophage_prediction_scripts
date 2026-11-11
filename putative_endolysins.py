from Bio import SeqIO
from collections import defaultdict
import os
import csv


gene_counts = defaultdict(int)

with open('cecorum_endolysins.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['Endolysin Name', 'Prophage Name', 'BiosampleID', 'Gene Before', 'Gene After',
                         'Endolysin Sequence'])
with open('putative_cecorum_endolysins.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['Gene Name', 'Prophage Name', 'BiosampleID', 'PHROGs Annotation', 'Gene Number', 'Holin Gene Number',
                         'Protein Sequence'])
with open('all_cecorum_endolysins.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['Gene Name', 'Putative or Annotated', 'Other Name','Protein Sequence'])

os.chdir('cecorum_prophages')
for file in os.listdir():
    if file.endswith('.gbk'):
        prophage_name = file.split('.')[0]
        biosampleid = prophage_name.split('_')[0]
        endolysin_count = 0
        gene_count = 0
        putative_endolysin_count = 0
        all_endolysin_count = 0
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
                        all_endolysin_count += 1
                        try:
                            next_gene = seq_record.features[gene_count + 1].qualifiers["product"][0]
                        except:
                            continue
                        endolysin_sequence = seq_record.features[gene_count].qualifiers["translation"][0]
                        with open('../cecorum_endolysins.tsv', 'a') as outfile:
                            tsv_writer = csv.writer(outfile, delimiter='\t')
                            tsv_writer.writerow(
                                [f'{prophage_name}_e{endolysin_count}', f'{prophage_name}', f'{biosampleid}',
                                 f'{previous_gene}', f'{next_gene}', f'{endolysin_sequence}'])
                        with open('../all_cecorum_endolysins.tsv', 'a') as outfile:
                            tsv_writer = csv.writer(outfile, delimiter='\t')
                            tsv_writer.writerow(
                                [f'{prophage_name}_{all_endolysin_count}', 'Annotated', f'{prophage_name}_e{endolysin_count}',f'{endolysin_sequence}'])
                        with open('../cecorum_endolysins.faa', 'a') as outfile:
                            outfile.write(f'>{prophage_name}_{all_endolysin_count}' + '\n' + endolysin_sequence + '\n')
                    elif current_gene == "holin":
                        holin_gene_number = gene_count
                        try:
                            two_genes_before = seq_record.features[gene_count - 2].qualifiers["product"][0]
                        except:
                            continue
                        try:
                            one_gene_before = seq_record.features[gene_count - 1].qualifiers["product"][0]
                        except:
                            continue
                        try:
                            one_gene_after = seq_record.features[gene_count + 1].qualifiers["product"][0]
                        except:
                            continue
                        try:
                            two_genes_after = seq_record.features[gene_count + 2].qualifiers["product"][0]
                        except:
                            continue
                        if not two_genes_before == "endolysin" and not two_genes_after == "endolysin" and not one_gene_before == "endolysin" and not one_gene_after == "endolysin" and not two_genes_after == "endolysin":
                            if two_genes_before == "hypothetical protein":
                                putative_endolysin_sequence = seq_record.features[gene_count - 2].qualifiers["translation"][0]
                                putative_endolysin_count += 1
                                all_endolysin_count += 1
                                with open('../putative_cecorum_endolysins.tsv', 'a') as outfile:
                                    tsv_writer = csv.writer(outfile, delimiter='\t')
                                    tsv_writer.writerow(
                                        [f'{prophage_name}_pe{putative_endolysin_count}', f'{prophage_name}', f'{biosampleid}',
                                         f'{two_genes_before}', gene_count - 2, f'{gene_count}', f'{putative_endolysin_sequence}'])
                                with open('../all_cecorum_endolysins.tsv', 'a') as outfile:
                                    tsv_writer = csv.writer(outfile, delimiter='\t')
                                    tsv_writer.writerow(
                                        [f'{prophage_name}_{all_endolysin_count}', 'Putative', f'{prophage_name}_pe{putative_endolysin_count}',
                                         f'{putative_endolysin_sequence}'])
                                with open('../cecorum_endolysins.faa', 'a') as outfile:
                                    outfile.write(f'>{prophage_name}_{all_endolysin_count}' + '\n' + putative_endolysin_sequence + '\n')
                            if one_gene_before == "hypothetical protein":
                                putative_endolysin_sequence = seq_record.features[gene_count - 1].qualifiers["translation"][0]
                                putative_endolysin_count += 1
                                all_endolysin_count += 1
                                with open('../putative_cecorum_endolysins.tsv', 'a') as outfile:
                                    tsv_writer = csv.writer(outfile, delimiter='\t')
                                    tsv_writer.writerow(
                                        [f'{prophage_name}_pe{putative_endolysin_count}', f'{prophage_name}', f'{biosampleid}',
                                         f'{one_gene_before}', gene_count - 1, f'{gene_count}', f'{putative_endolysin_sequence}'])
                                with open('../all_cecorum_endolysins.tsv', 'a') as outfile:
                                    tsv_writer = csv.writer(outfile, delimiter='\t')
                                    tsv_writer.writerow(
                                        [f'{prophage_name}_{all_endolysin_count}', 'Putative',f'{prophage_name}_pe{putative_endolysin_count}',
                                         f'{putative_endolysin_sequence}'])
                                with open('../cecorum_endolysins.faa', 'a') as outfile:
                                    outfile.write(
                                        f'>{prophage_name}_{all_endolysin_count}' + '\n' + putative_endolysin_sequence + '\n')
                            if one_gene_after == "hypothetical protein":
                                putative_endolysin_sequence = seq_record.features[gene_count + 1].qualifiers["translation"][0]
                                putative_endolysin_count += 1
                                all_endolysin_count += 1
                                with open('../putative_cecorum_endolysins.tsv', 'a') as outfile:
                                    tsv_writer = csv.writer(outfile, delimiter='\t')
                                    tsv_writer.writerow(
                                        [f'{prophage_name}_pe{putative_endolysin_count}', f'{prophage_name}', f'{biosampleid}',
                                         f'{one_gene_after}', gene_count + 1, f'{gene_count}', f'{putative_endolysin_sequence}'])
                                with open('../all_cecorum_endolysins.tsv', 'a') as outfile:
                                    tsv_writer = csv.writer(outfile, delimiter='\t')
                                    tsv_writer.writerow(
                                        [f'{prophage_name}_{all_endolysin_count}', 'Putative',f'{prophage_name}_pe{putative_endolysin_count}',
                                         f'{putative_endolysin_sequence}'])
                                with open('../cecorum_endolysins.faa', 'a') as outfile:
                                    outfile.write(
                                        f'>{prophage_name}_{all_endolysin_count}' + '\n' + putative_endolysin_sequence + '\n')
                            if two_genes_after == "hypothetical protein":
                                putative_endolysin_sequence = seq_record.features[gene_count + 2].qualifiers["translation"][0]
                                putative_endolysin_count += 1
                                all_endolysin_count += 1
                                with open('../putative_cecorum_endolysins.tsv', 'a') as outfile:
                                    tsv_writer = csv.writer(outfile, delimiter='\t')
                                    tsv_writer.writerow(
                                        [f'{prophage_name}_pe{putative_endolysin_count}', f'{prophage_name}', f'{biosampleid}',
                                         f'{two_genes_after}', gene_count + 2, f'{gene_count}', f'{putative_endolysin_sequence}'])
                                with open('../all_cecorum_endolysins.tsv', 'a') as outfile:
                                    tsv_writer = csv.writer(outfile, delimiter='\t')
                                    tsv_writer.writerow(
                                        [f'{prophage_name}_{all_endolysin_count}', 'Putative',f'{prophage_name}_pe{putative_endolysin_count}',
                                         f'{putative_endolysin_sequence}'])
                                with open('../cecorum_endolysins.faa', 'a') as outfile:
                                    outfile.write(
                                        f'>{prophage_name}_{all_endolysin_count}' + '\n' + putative_endolysin_sequence + '\n')
                gene_count += 1

