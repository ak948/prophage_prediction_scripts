"""
Output file with rows for each prophage name
Take BiosampleID in results.tsv and prophage number from it, unzip that one, and read phageboost_results/prophages.txt
For prophage in prophage_dataframe.tsv (until prophagenumber > prophage number from results),
read 2nd Name until it finds final_{biosampleid}.fna, for that row take the value in Comment (Split between 2nd and 3rd _ for Genomad or VIBRANT)
From value in Comment, take _prophage{number} from final{biosampleid}_{software}_prophage{number}.fna, then find _prophage{number} in Discard 1, and save comment in Comment
Read 1st Change until find _prophage{number}, and save value in Original Filename and take _prophage{number} from it as {originalnumber}

Read phages_{BiosampleID}.gff, ignore first row, from attributes find phage_id=phage{originalnumber}, from that row take value in #seqid, start, end
Read bacteria genome and find contig matching #seqid, find start base and take 500 bases before, find end base and take 500 bases after

Output row should have Prophage Name, Genomad or VIBRANT, Why other one discarded, Renamed Phageboost Name, Original Phageboost Name, Phageboost Start,
Phageboost End, Phageboost Seqid, Nucleotides Before, Nucleotides After
"""

import os
import pandas as pd
import subprocess
import csv
from Bio import SeqIO

with open('prophage_location_info.tsv', 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['Prophage Name', 'Software Used', 'Why Other One Discarded', 'Renamed Phageboost Name',
                         'Original Phageboost Name', 'Phageboost Start', 'Phageboost End', 'Phageboost SeqID',
                         'Nucleotides Before', 'Nucleotides After', 'Prophage Start Location', 'Prophage End Location'])

df = pd.read_csv('results.tsv', sep='\t')
genomeid = df['BiosampleID']


row_num = 0
for biosampleid in genomeid:
    prophage_count = 1
    prophage_number = df['Number of Prophages'][row_num]
    os.chdir('Completed')
    subprocess.run(f'tar -zxvf {biosampleid}.tar.gz {biosampleid}', shell=True)
    os.chdir('..')
    while prophage_count <= prophage_number:
        if os.path.isfile(f'Completed/{biosampleid}/phageboost_results/prophages.txt') and os.path.isfile(f'Completed/{biosampleid}/phageboost_results/phages_{biosampleid}.gff'):
            with open(f'Completed/{biosampleid}/phageboost_results/prophages.txt', 'r') as file:
                reader = csv.DictReader(file, delimiter='\t')
                for row in reader:
                    if row['2nd Name'] == f'final_{biosampleid}_prophage{prophage_count}.fna':
                        second_name = row['Comment']
                        software_used = second_name.split('_')[2]
                        second_half = second_name.split('_')[3]
                        second_prophage_number = second_half.split('.')[0]
                        break
                    elif row['2nd Name'] == f'final_{biosampleid}_p{prophage_count}.fna':
                        second_name = row['Comment']
                        software_used = second_name.split('_')[2]
                        second_half = second_name.split('_')[3]
                        second_prophage_number = second_half.split('.')[0]
                        break
                    else:
                        continue
            with open(f'Completed/{biosampleid}/phageboost_results/prophages.txt', 'r') as file:
                reader2 = csv.DictReader(file, delimiter='\t')
                for row2 in reader2:
                    if row2['Comment'] == 'Both Empty':
                        if row2['Discard 2'] == f'final_{biosampleid}_genomad_{second_prophage_number}.fna':
                            comment = row2['Comment']
                            break
                    else:
                        if row2['Keep'] == f'final_{biosampleid}_{software_used}_{second_prophage_number}.fna':
                            comment = row2['Discard 2']
                            break

            with open(f'Completed/{biosampleid}/phageboost_results/prophages.txt', 'r') as file:
                reader3 = csv.DictReader(file, delimiter='\t')
                for row3 in reader3:
                    if row3['1st Change'] == f'final_{biosampleid}_phageboost_{second_prophage_number}.fasta':
                        renamed_name = row3['1st Change']
                        original_name = row3['Original Filename']
                        first_split = original_name.split('_')[1]
                        original_prophage_number = first_split.split('.')[0]
                        break
            with open(f'Completed/{biosampleid}/phageboost_results/phages_{biosampleid}.gff', 'r') as file:
                # Skip metadata lines starting with '#'
                lines = [line for line in file if not line.startswith('##')]
                # Create a DictReader using the actual header from the GFF (2nd # line)
                reader = csv.DictReader(lines, delimiter='\t')
                for row in reader:
                    attributes = row['attributes']
                    phageid = attributes.split('=')[2]
                    if phageid == original_prophage_number:
                        contig_id = row['#seqid']
                        prophage_start = int(row['start'])
                        prophage_end = int(row['end'])
                        break

            sequence_before = None
            sequence_after = None
            current_genome_length = 0
            for seq_record in SeqIO.parse(f'BacterialGenomes/{biosampleid}.fna', 'fasta'):
                if seq_record.id == contig_id:
                    seq_len = len(seq_record.seq)
                    flank_length = 2000
                    start_index = max(prophage_start - flank_length - 1, 0)
                    end_index = prophage_start - 1
                    sequence_before = seq_record.seq[start_index:end_index]
                    start_index = prophage_end  # because GFF is inclusive, and Python indexing starts at 0
                    end_index = min(prophage_end + flank_length, seq_len)  # prevent out-of-bounds

                    sequence_after = seq_record.seq[start_index:end_index]
                    prophage_start_location = current_genome_length + prophage_start
                    prophage_end_location = current_genome_length + prophage_end
                    break
                else:
                    current_genome_length += len(seq_record.seq)

            with open('prophage_location_info.tsv', 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow(
                    [f'{biosampleid}_p{prophage_count}', f'{software_used}', f'{comment}', f'{renamed_name}', f'{original_name}', f'{prophage_start}',
                     f'{prophage_end}', f'{contig_id}', f'{sequence_before}', f'{sequence_after}',
                     f'{prophage_start_location}', f'{prophage_end_location}'])

            prophage_count += 1
        else:
            with open('prophage_location_info.tsv', 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow(
                    [f'{biosampleid}_p{prophage_count}', 'No Info', f'No Info', f'No Info', f'No Info', f'No Info',
                     f'No Info', f'No Info', f'No Info', f'No Info', f'No Info', f'No Info'])
            prophage_count += 1
    row_num +=1
    os.chdir('Completed')
    subprocess.run(f'rm -r {biosampleid}', shell=True)
    os.chdir('..')
