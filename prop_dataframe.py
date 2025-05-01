

with open('prophage_dataframe.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Prophage Name', 'BiosampleID', 'Genome Length', 'Coding Sequences', 'GC%',
                         'Virulence Gene Present',
                         'Virulence Gene Count', "Virulence Gene 'ecbA/fss3' Count", "Virulence Gene 'srtC' Count",
                         "Virulence Gene 'ebpC' Count", "Virulence Gene 'ebpB' Count", "Virulence Gene 'ebpA' Count",
                         "Virulence Gene 'ace' Count", "Virulence Gene 'prgB/asc10' Count",
                         "Virulence Gene 'asa1' Count", "Virulence Gene 'EF0485' Count", "Virulence Gene 'acm' Count",
                         "Virulence Gene 'cylR2' Count", "Virulence Gene 'cpsI' Count", "Virulence Gene 'cylS' Count",
                         "Virulence Gene 'cpsC' Count", "Virulence Gene 'cpsD' Count", "Virulence Gene 'cpsE' Count",
                         "Virulence Gene 'cpsG' Count", "Virulence Gene 'cpsH' Count", "Virulence Gene 'cylR1' Count",
                         "Virulence Gene 'cpsJ' Count", "Virulence Gene 'cpsK' Count", "Virulence Gene 'cylL' Count",
                         "Virulence Gene 'cylB' Count", "Virulence Gene 'fsrB' Count", "Virulence Gene 'fsrA' Count",
                         "Virulence Gene 'cylM' Count", "Virulence Gene 'cylA' Count", "Virulence Gene 'EF3023' Count",
                         "Virulence Gene 'cpsF' Count", "Virulence Gene 'efaA' Count", 'ARG Present', 'ARG Count',
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
                         "ARG Family 'vanY; glycopeptide resistance gene cluster' Count", 'Cluster', 'Phage Genus',
                         'Phage Species', 'Host Genbank Species', 'Sample Collection Date', 'Sample Collection Year',
                         'Sample Collection Location', 'Sample Collection Country', 'Sample Isolation Source',
                         'Number of Prophages in Host', 'Average Prophage Length in Host',
                         'Minimum Prophage Length in Host', 'Maximum Prophage Length in Host',
                         'Average number of coding sequences in prophages in Host',
                         'Minimum number of coding sequences in prophages in host',
                         'Maximum number of coding sequences in prophages in host',
                         'Average Prophage GC% in Host', 'Minimum Prophage GC% in Host', 'Maximum Prophage GC% in Host',
                         'Number of prophages containing virulence genes in host',
                         'Total number of virulence genes in prophages in host',
                         'Number of prophages containing AMR Genes in host',
                         'Total number of AMR Genes in prophages in host'])

df = pd.read_csv('prophage_results.tsv', sep='\t')
prophage_name = df['Prophage Name']

for name in prophage_name:
    vir_gene_1 = 0
    vir_gene_2 = 0
    vir_gene_3 = 0
    vir_gene_4 = 0
    vir_gene_5 = 0
    vir_gene_6 = 0
    vir_gene_7 = 0
    vir_gene_8 = 0
    vir_gene_9 = 0
    vir_gene_10 = 0
    vir_gene_11 = 0
    vir_gene_12 = 0
    vir_gene_13 = 0
    vir_gene_14 = 0
    vir_gene_15 = 0
    vir_gene_16 = 0
    vir_gene_17 = 0
    vir_gene_18 = 0
    vir_gene_19 = 0
    vir_gene_20 = 0
    vir_gene_21 = 0
    vir_gene_22 = 0
    vir_gene_23 = 0
    vir_gene_24 = 0
    vir_gene_25 = 0
    vir_gene_26 = 0
    vir_gene_27 = 0
    vir_gene_28 = 0
    vir_gene_29 = 0
    vir_gene_30 = 0
    amr_gene_1 = 0
    amr_gene_2 = 0
    amr_gene_3 = 0
    amr_gene_4 = 0
    amr_gene_5 = 0
    amr_gene_6 = 0
    amr_gene_7 = 0
    amr_gene_8 = 0
    amr_gene_9 = 0
    amr_gene_10 = 0
    amr_gene_11 = 0
    amr_gene_12 = 0
    amr_gene_13 = 0
    amr_gene_14 = 0
    amr_gene_15 = 0
    amr_gene_16 = 0
    amr_gene_17 = 0
    amr_gene_18 = 0
    amr_gene_19 = 0
    amr_gene_20 = 0
    amr_gene_21 = 0
    amr_gene_22 = 0
    amr_gene_23 = 0
    amr_gene_24 = 0
    amr_gene_25 = 0
    amr_gene_26 = 0

    with open('prophage_results.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['Prophage Name'] == name:
                biosampleid = row['Biosampleid']
                genome_length = row['Genome Length']
                coding_sequences = row['Coding Sequences']
                gcpercent = row['GC%']
                virulence_gene_count = row['Virulence Gene Count']
                if int(virulence_gene_count) > 0:
                    virulence_present = '1'
                else:
                    virulence_present = '0'
                amr_count = row['AMR Count']
                if int(amr_count) > 0:
                    amr_present = '1'
                else:
                    amr_present = '0'
    with open('virulence_genes.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['Sequence'] == name:
                if row['Gene'] == 'ecbA/fss3':
                    vir_gene_1 += 1
                elif row['Gene'] == 'srtC':
                    vir_gene_2 += 1
                elif row['Gene'] == 'ebpC':
                    vir_gene_3 += 1
                elif row['Gene'] == 'ebpB':
                    vir_gene_4 += 1
                elif row['Gene'] == 'ebpA':
                    vir_gene_5 += 1
                elif row['Gene'] == 'ace':
                    vir_gene_6 += 1
                elif row['Gene'] == 'prgB/asc10':
                    vir_gene_7 += 1
                elif row['Gene'] == 'asa1':
                    vir_gene_8 += 1
                elif row['Gene'] == 'EF0485':
                    vir_gene_9 += 1
                elif row['Gene'] == 'acm':
                    vir_gene_10 += 1
                elif row['Gene'] == 'cylR2':
                    vir_gene_11 += 1
                elif row['Gene'] == 'cpsI':
                    vir_gene_12 += 1
                elif row['Gene'] == 'cylS':
                    vir_gene_13 += 1
                elif row['Gene'] == 'cpsC':
                    vir_gene_14 += 1
                elif row['Gene'] == 'cpsD':
                    vir_gene_15 += 1
                elif row['Gene'] == 'cpsE':
                    vir_gene_16 += 1
                elif row['Gene'] == 'cpsG':
                    vir_gene_17 += 1
                elif row['Gene'] == 'cpsH':
                    vir_gene_18 += 1
                elif row['Gene'] == 'cylR1':
                    vir_gene_19 += 1
                elif row['Gene'] == 'cpsJ':
                    vir_gene_20 += 1
                elif row['Gene'] == 'cpsK':
                    vir_gene_21 += 1
                elif row['Gene'] == 'cylL':
                    vir_gene_22 += 1
                elif row['Gene'] == 'cylB':
                    vir_gene_23 += 1
                elif row['Gene'] == 'fsrB':
                    vir_gene_24 += 1
                elif row['Gene'] == 'fsrA':
                    vir_gene_25 += 1
                elif row['Gene'] == 'cylM':
                    vir_gene_26 += 1
                elif row['Gene'] == 'cylA':
                    vir_gene_27 += 1

                elif row['Gene'] == 'EF3023':
                    vir_gene_28 += 1
                elif row['Gene'] == 'cpsF':
                    vir_gene_29 += 1
                elif row['Gene'] == 'efaA':
                    vir_gene_30 += 1
    with open('amr_genes.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{name}_' in row['Source Sequence']:
                if row['AMR Gene Family'] == 'ADC beta-lactamases pending classification for carbapenemase activity':
                    amr_gene_1 += 1
                elif row['AMR Gene Family'] == 'aminoglycoside bifunctional resistance protein':
                    amr_gene_2 += 1
                elif row['AMR Gene Family'] == 'ANT(6)':
                    amr_gene_3 += 1
                elif row['AMR Gene Family'] == 'ANT(9)':
                    amr_gene_4 += 1
                elif row['AMR Gene Family'] == "APH(3')":
                    amr_gene_5 += 1
                elif row['AMR Gene Family'] == 'ATP-binding cassette (ABC) antibiotic efflux pump':
                    amr_gene_6 += 1
                elif row['AMR Gene Family'] == 'Cfr 23S ribosomal RNA methyltransferase':
                    amr_gene_7 += 1
                elif row['AMR Gene Family'] == 'Erm 23S ribosomal RNA methyltransferase':
                    amr_gene_8 += 1
                elif row['AMR Gene Family'] == 'fosfomycin thiol transferase':
                    amr_gene_9 += 1
                elif row['AMR Gene Family'] == 'glycopeptide resistance gene cluster; vanR':
                    amr_gene_10 += 1
                elif row['AMR Gene Family'] == 'lincosamide nucleotidyltransferase (LNU)':
                    amr_gene_11 += 1
                elif row['AMR Gene Family'] == 'lsa-type ABC-F protein':
                    amr_gene_12 += 1
                elif row['AMR Gene Family'] == 'major facilitator superfamily (MFS) antibiotic efflux pump':
                    amr_gene_13 += 1
                elif row['AMR Gene Family'] == 'Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal protection proteins':
                    amr_gene_14 += 1
                elif row['AMR Gene Family'] == 'msr-type ABC-F protein':
                    amr_gene_15 += 1
                elif row['AMR Gene Family'] == 'multidrug and toxic compound extrusion (MATE) transporter':
                    amr_gene_16 += 1
                elif row['AMR Gene Family'] == 'small multidrug resistance (SMR) antibiotic efflux pump':
                    amr_gene_17 += 1
                elif row['AMR Gene Family'] == 'streptothricin acetyltransferase (SAT)':
                    amr_gene_18 += 1
                elif row['AMR Gene Family'] == 'TEM beta-lactamase':
                    amr_gene_19 += 1
                elif row['AMR Gene Family'] == 'tetracycline-resistant ribosomal protection protein':
                    amr_gene_20 += 1
                elif row['AMR Gene Family'] == 'trimethoprim resistant dihydrofolate reductase dfr':
                    amr_gene_21 += 1
                elif row['AMR Gene Family'] == 'vanH; glycopeptide resistance gene cluster':
                    amr_gene_22 += 1
                elif row['AMR Gene Family'] == 'vanS; glycopeptide resistance gene cluster':
                    amr_gene_23 += 1
                elif row['AMR Gene Family'] == 'vanW; glycopeptide resistance gene cluster':
                    amr_gene_24 += 1
                elif row['AMR Gene Family'] == 'vanX; glycopeptide resistance gene cluster':
                    amr_gene_25 += 1
                elif row['AMR Gene Family'] == 'vanY; glycopeptide resistance gene cluster':
                    amr_gene_26 += 1
    with open('my_clusters.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t', fieldnames=['column1', 'column2'])
        for row in reader:
            if f'{name},' in row['column2']:
                cluster_rep = row['column1']
            elif f'{name}"' in row['column2']:
                cluster_rep = row['column1']
            elif row['column2'] == name:
                cluster_rep = row['column1']
    with open('Summary_taxonomy.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{cluster_rep}' in row['Genome']:
                cluster = row['Cluster']
                phage_genus = row['Genus']
                phage_species = row['Species']
    with open('results.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if f'{biosampleid}' in row['Biosampleid']:
                genbank_species = row['Genbank Info Species']
                collection_date = row['Collection Date']
                collection_year = row['Collection Year']
                collection_location = row['Collection Location']
                collection_country = row['Collection Country']
                isolation_source = row['Isolation Source']
                host_prophages = row['Number of Prophages']
                host_average_prophage_length = row['Average Prophage Length']
                host_minimum_prophage_length = row['Minimum Prophage Length']
                host_maximum_prophage_length = row['Maximum Prophage Length']
                host_average_prophage_cds = row['Average number of prophage coding sequences']
                host_minimum_prophage_cds = row['Minimum number of prophage coding sequences']
                host_maximum_prophage_cds = row['Maximum number of prophage coding sequences']
                host_average_prophage_gc = row['Average Prophage GC%']
                host_minimum_prophage_gc = row['Minimum Prophage GC%']
                host_maximum_prophage_gc = row['Maximum Prophage GC%']
                host_prophage_vir = row['Number of prophages containing virulence genes']
                host_prophage_vir_no = row['Total number of virulence genes in prophages']
                host_prophage_amr = row['Number of prophages containing AMR Genes']
                host_prophage_amr_no = row['Total number of AMR Genes in prophages']
    with open('prophage_dataframe.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow([f'{name}', f'{biosampleid}', f'{genome_length}', f'{coding_sequences}',
                             f'{gcpercent}', f'{virulence_present}', f'{virulence_gene_count}', f'{vir_gene_1}',
                             f'{vir_gene_2}', f'{vir_gene_3}', f'{vir_gene_4}', f'{vir_gene_5}', f'{vir_gene_6}',
                             f'{vir_gene_7}', f'{vir_gene_8}', f'{vir_gene_9}', f'{vir_gene_10}', f'{vir_gene_11}',
                             f'{vir_gene_12}', f'{vir_gene_13}', f'{vir_gene_14}', f'{vir_gene_15}', f'{vir_gene_16}',
                             f'{vir_gene_17}', f'{vir_gene_18}', f'{vir_gene_19}', f'{vir_gene_20}', f'{vir_gene_21}',
                             f'{vir_gene_22}', f'{vir_gene_23}', f'{vir_gene_24}', f'{vir_gene_25}', f'{vir_gene_26}',
                             f'{vir_gene_27}', f'{vir_gene_28}', f'{vir_gene_29}', f'{vir_gene_30}', f'{amr_present}',
                             f'{amr_count}', f'{amr_gene_1}', f'{amr_gene_2}', f'{amr_gene_3}', f'{amr_gene_4}',
                             f'{amr_gene_5}', f'{amr_gene_6}', f'{amr_gene_7}', f'{amr_gene_8}', f'{amr_gene_9}',
                             f'{amr_gene_10}', f'{amr_gene_11}', f'{amr_gene_12}', f'{amr_gene_13}', f'{amr_gene_14}',
                             f'{amr_gene_15}', f'{amr_gene_16}', f'{amr_gene_17}', f'{amr_gene_18}', f'{amr_gene_19}',
                             f'{amr_gene_20}', f'{amr_gene_21}', f'{amr_gene_22}', f'{amr_gene_23}', f'{amr_gene_24}',
                             f'{amr_gene_25}', f'{amr_gene_26}', f'{cluster}', f'{phage_genus}', f'{phage_species}',
                             f'{genbank_species}', f'{collection_date}', f'{collection_year}', f'{collection_location}',
                             f'{collection_country}', f'{isolation_source}', f'{host_prophages}',
                             f'{host_average_prophage_length}', f'{host_minimum_prophage_length}',
                             f'{host_maximum_prophage_length}', f'{host_average_prophage_cds}',
                             f'{host_minimum_prophage_cds}', f'{host_maximum_prophage_cds}',
                             f'{host_average_prophage_gc}', f'{host_minimum_prophage_gc}',
                             f'{host_maximum_prophage_gc}', f'{host_prophage_vir}', f'{host_prophage_vir_no}',
                             f'{host_prophage_amr}', f'{host_prophage_amr_no}'])
