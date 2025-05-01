import os
import pandas as pd
import csv

with open('removedgenomes.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(
        ['BiosampleID',	'Genbank Info Species',	'Collection Date',	'Collection Year',	'Collection Location',
         'Collection Country',	'Isolation Source',	'Number of Prophages',	'Average Prophage Length',
         'Minimum Prophage Length',	'Maximum Prophage Length',	'Average number of prophage coding sequences',
         'Minimum number of prophage coding sequences',	'Maximum number of prophage coding sequences',
         'Average Prophage GC%',	'Minimum Prophage GC%',	'Maximum Prophage GC%',
         'Number of prophages containing virulence genes',	'Total number of virulence genes in prophages',
         'Number of prophages containing AMR Genes',	'Total number of AMR Genes in prophages',	'Classification',
         'Checkm Marker lineage',	'Checkm # genomes',	'Checkm # markers',	'Checkm # marker sets',
         'Checkm 0 Single Copies',	'Checkm 1 Single Copies',	'Checkm 2 Single Copies',	'Checkm 3 Single Copies',
         'Checkm 4 Single Copies', 'Checkm 5+ Single Copies',	'Checkm Completeness',	'Checkm Contamination',
         'Checkm Strain heterogeneity'])

with open('removedprophages.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(["Prophage Name", "BiosampleID", "Genome Length", "Coding Sequences", "GC%",
                         "Virulence Gene Present", "Virulence Gene Count", "Virulence Gene 'ecbA/fss3' Count",
                         "Virulence Gene 'srtC' Count", "Virulence Gene 'ebpC' Count", "Virulence Gene 'ebpB' Count",
                         "Virulence Gene 'ebpA' Count", "Virulence Gene 'ace' Count",
                         "Virulence Gene 'prgB/asc10' Count", "Virulence Gene 'asa1' Count",
                         "Virulence Gene 'EF0485' Count", "Virulence Gene 'acm' Count", "Virulence Gene 'cylR2' Count",
                         "Virulence Gene 'cpsI' Count", "Virulence Gene 'cylS' Count", "Virulence Gene 'cpsC' Count",
                         "Virulence Gene 'cpsD' Count", "Virulence Gene 'cpsE' Count", "Virulence Gene 'cpsG' Count",
                         "Virulence Gene 'cpsH' Count", "Virulence Gene 'cylR1' Count", "Virulence Gene 'cpsJ' Count",
                         "Virulence Gene 'cpsK' Count", "Virulence Gene 'cylL' Count", "Virulence Gene 'cylB' Count",
                         "Virulence Gene 'fsrB' Count", "Virulence Gene 'fsrA' Count", "Virulence Gene 'cylM' Count",
                         "Virulence Gene 'cylA' Count", "Virulence Gene 'EF3023' Count", "Virulence Gene 'cpsF' Count",
                         "Virulence Gene 'efaA' Count", "ARG Present", "ARG Count",
                         "ARG Family 'ADC beta-lactamases pending classification for carbapenemase activity' Count",
                         "ARG Family 'aminoglycoside bifunctional resistance protein' Count",
                         "ARG Family 'ANT(6)' Count", "ARG Family 'ANT(9)' Count", "ARG Family 'APH(3')' Count",
                         "ARG Family 'ATP-binding cassette (ABC) antibiotic efflux pump' Count",
                         "ARG Family 'Cfr 23S ribosomal RNA methyltransferase' Count",
                         "ARG Family 'Erm 23S ribosomal RNA methyltransferase' Count",
                         "ARG Family 'fosfomycin thiol transferase' Count",
                         "ARG Family 'glycopeptide resistance gene cluster; vanR' Count",
                         "ARG Family 'lincosamide nucleotidyltransferase (LNU)' Count",
                         "ARG Family 'lsa-type ABC-F protein' Count",
                         "ARG Family 'major facilitator superfamily (MFS) antibiotic efflux pump' Count",
                         "ARG Family 'Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal protection proteins' Count",
                         "ARG Family 'msr-type ABC-F protein' Count",
                         "ARG Family 'multidrug and toxic compound extrusion (MATE) transporter' Count",
                         "ARG Family 'small multidrug resistance (SMR) antibiotic efflux pump' Count",
                         "ARG Family 'streptothricin acetyltransferase (SAT)' Count",
                         "ARG Family 'TEM beta-lactamase' Count",
                         "ARG Family 'tetracycline-resistant ribosomal protection protein' Count",
                         "ARG Family 'trimethoprim resistant dihydrofolate reductase dfr' Count",
                         "ARG Family 'vanH; glycopeptide resistance gene cluster' Count",
                         "ARG Family 'vanS; glycopeptide resistance gene cluster' Count",
                         "ARG Family 'vanW; glycopeptide resistance gene cluster' Count",
                         "ARG Family 'vanX; glycopeptide resistance gene cluster' Count",
                         "ARG Family 'vanY; glycopeptide resistance gene cluster' Count", "Cluster", "Phage Genus",
                         "Phage Species", "Host Genbank Species", "Sample Collection Date", "Sample Collection Year",
                         "Sample Collection Location", "Sample Collection Country", "Sample Isolation Source",
                         "Number of Prophages in Host", "Average Prophage Length in Host",
                         "Minimum Prophage Length in Host", "Maximum Prophage Length in Host",
                         "Average number of coding sequences in prophages in Host",
                         "Minimum number of coding sequences in prophages in host",
                         "Maximum number of coding sequences in prophages in host", "Average Prophage GC% in Host",
                         "Minimum Prophage GC% in Host", "Maximum Prophage GC% in Host",
                         "Number of prophages containing virulence genes in host",
                         "Total number of virulence genes in prophages in host",
                         "Number of prophages containing AMR Genes in host",
                         "Total number of AMR Genes in prophages in host", "Cluster Representative", "Cluster Quality",
                         "VC", "Genomes in VC", "Representative Genomes in VC", "Total Prophage Genomes in VC",
                         "Classification", "Checkm Marker lineage", "Checkm # genomes", "Checkm # markers",
                         "Checkm # marker sets", "Checkm 0 Single Copies", "Checkm 1 Single Copies",
                         "Checkm 2 Single Copies", "Checkm 3 Single Copies", "Checkm 4 Single Copies",
                         "Checkm 5+ Single Copies", "Checkm Completeness", "Checkm Contamination",
                         "Checkm Strain heterogeneity"])

with open('filteredprophages.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(["Prophage Name", "BiosampleID", "Genome Length", "Coding Sequences", "GC%",
                         "Virulence Gene Present", "Virulence Gene Count", "Virulence Gene 'ecbA/fss3' Count",
                         "Virulence Gene 'srtC' Count", "Virulence Gene 'ebpC' Count", "Virulence Gene 'ebpB' Count",
                         "Virulence Gene 'ebpA' Count", "Virulence Gene 'ace' Count",
                         "Virulence Gene 'prgB/asc10' Count", "Virulence Gene 'asa1' Count",
                         "Virulence Gene 'EF0485' Count", "Virulence Gene 'acm' Count", "Virulence Gene 'cylR2' Count",
                         "Virulence Gene 'cpsI' Count", "Virulence Gene 'cylS' Count", "Virulence Gene 'cpsC' Count",
                         "Virulence Gene 'cpsD' Count", "Virulence Gene 'cpsE' Count", "Virulence Gene 'cpsG' Count",
                         "Virulence Gene 'cpsH' Count", "Virulence Gene 'cylR1' Count", "Virulence Gene 'cpsJ' Count",
                         "Virulence Gene 'cpsK' Count", "Virulence Gene 'cylL' Count", "Virulence Gene 'cylB' Count",
                         "Virulence Gene 'fsrB' Count", "Virulence Gene 'fsrA' Count", "Virulence Gene 'cylM' Count",
                         "Virulence Gene 'cylA' Count", "Virulence Gene 'EF3023' Count", "Virulence Gene 'cpsF' Count",
                         "Virulence Gene 'efaA' Count", "ARG Present", "ARG Count",
                         "ARG Family 'ADC beta-lactamases pending classification for carbapenemase activity' Count",
                         "ARG Family 'aminoglycoside bifunctional resistance protein' Count",
                         "ARG Family 'ANT(6)' Count", "ARG Family 'ANT(9)' Count", "ARG Family 'APH(3')' Count",
                         "ARG Family 'ATP-binding cassette (ABC) antibiotic efflux pump' Count",
                         "ARG Family 'Cfr 23S ribosomal RNA methyltransferase' Count",
                         "ARG Family 'Erm 23S ribosomal RNA methyltransferase' Count",
                         "ARG Family 'fosfomycin thiol transferase' Count",
                         "ARG Family 'glycopeptide resistance gene cluster; vanR' Count",
                         "ARG Family 'lincosamide nucleotidyltransferase (LNU)' Count",
                         "ARG Family 'lsa-type ABC-F protein' Count",
                         "ARG Family 'major facilitator superfamily (MFS) antibiotic efflux pump' Count",
                         "ARG Family 'Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal protection proteins' Count",
                         "ARG Family 'msr-type ABC-F protein' Count",
                         "ARG Family 'multidrug and toxic compound extrusion (MATE) transporter' Count",
                         "ARG Family 'small multidrug resistance (SMR) antibiotic efflux pump' Count",
                         "ARG Family 'streptothricin acetyltransferase (SAT)' Count",
                         "ARG Family 'TEM beta-lactamase' Count",
                         "ARG Family 'tetracycline-resistant ribosomal protection protein' Count",
                         "ARG Family 'trimethoprim resistant dihydrofolate reductase dfr' Count",
                         "ARG Family 'vanH; glycopeptide resistance gene cluster' Count",
                         "ARG Family 'vanS; glycopeptide resistance gene cluster' Count",
                         "ARG Family 'vanW; glycopeptide resistance gene cluster' Count",
                         "ARG Family 'vanX; glycopeptide resistance gene cluster' Count",
                         "ARG Family 'vanY; glycopeptide resistance gene cluster' Count", "Cluster", "Phage Genus",
                         "Phage Species", "Host Genbank Species", "Sample Collection Date", "Sample Collection Year",
                         "Sample Collection Location", "Sample Collection Country", "Sample Isolation Source",
                         "Number of Prophages in Host", "Average Prophage Length in Host",
                         "Minimum Prophage Length in Host", "Maximum Prophage Length in Host",
                         "Average number of coding sequences in prophages in Host",
                         "Minimum number of coding sequences in prophages in host",
                         "Maximum number of coding sequences in prophages in host", "Average Prophage GC% in Host",
                         "Minimum Prophage GC% in Host", "Maximum Prophage GC% in Host",
                         "Number of prophages containing virulence genes in host",
                         "Total number of virulence genes in prophages in host",
                         "Number of prophages containing AMR Genes in host",
                         "Total number of AMR Genes in prophages in host", "Cluster Representative", "Cluster Quality",
                         "VC", "Genomes in VC", "Representative Genomes in VC", "Total Prophage Genomes in VC",
                         "Classification", "Checkm Marker lineage", "Checkm # genomes", "Checkm # markers",
                         "Checkm # marker sets", "Checkm 0 Single Copies", "Checkm 1 Single Copies",
                         "Checkm 2 Single Copies", "Checkm 3 Single Copies", "Checkm 4 Single Copies",
                         "Checkm 5+ Single Copies", "Checkm Completeness", "Checkm Contamination",
                         "Checkm Strain heterogeneity"])

id_to_search = "BiosampleID"
df = pd.read_csv('removed.txt', sep='\t')
biosample_ids = df[id_to_search]
biosampleid = df['BiosampleID']

new_row_num = 0
for id in biosampleid:
    row_num = 0
    with open('results.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if id == row['BiosampleID']:
                dg = pd.read_csv('removedgenomes.tsv', sep='\t', low_memory=False)
                dt = pd.read_csv('results.tsv', sep='\t', low_memory=False)
                dg.loc[new_row_num] = dt.loc[row_num]
                dg.to_csv('removedgenomes.tsv', sep='\t', index=False)
                dt.drop(row_num, inplace=True)
                new_row_num +=1
                dt.to_csv('results.tsv', sep='\t', index=False)
            row_num += 1

dt = pd.read_csv('prophage_dataframe.tsv', sep='\t', low_memory=False)
prophage_row_num = 0
while len(dt.index) > 0:
    if dt.loc[prophage_row_num][id_to_search] in biosample_ids.values:
        dg = pd.read_csv('removedprophages.tsv', sep='\t', low_memory=False)
        dg = dg.append(dt.loc[prophage_row_num])
        dg.to_csv('removedprophages.tsv', sep='\t', index=False)
        prophage_row_num +=1
    else:
        dh = pd.read_csv('filteredprophages.tsv', sep='\t', low_memory=False)
        dh = dh.append(dt.loc[prophage_row_num])
        dh.to_csv('filteredprophages.tsv', sep='\t', index=False)
        prophage_row_num += 1
    if prophage_row_num >= len(dt.index):
        break
