import pandas as pd
phage_df = pd.read_csv('14Apr2025_data.tsv', sep="\t", low_memory=False)
ictv_df = pd.read_csv('ICTV_Master_Species_List.tsv', sep="\t", low_memory=False)
taxonomy_levels = ['Genus', 'Subfamily', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom', 'Realm']
phage_df['Genome Composition'] = None

for level in taxonomy_levels:
    print(f"Trying to match on {level}...")
    if level in ictv_df.columns and level in phage_df.columns:
        mask = phage_df['Genome Composition'].isna()
        # Fix for duplicates in index
        temp_df = ictv_df[[level, 'Genome Composition']].dropna(subset=[level, 'Genome Composition'])

        # Option 1: drop duplicates by taxonomy level, keep first
        mapping = temp_df.drop_duplicates(subset=[level]).set_index(level)['Genome Composition']

        # Map genome composition from ICTV to phage_df
        phage_df.loc[mask, 'Genome Composition'] = phage_df.loc[mask, level].map(mapping)
    else:
        print(f"Column {level} missing in one of the files")
phage_df.to_excel('14Apr2025_data_with_composition.xlsx', index=False)
