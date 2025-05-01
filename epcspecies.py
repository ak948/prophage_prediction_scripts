# Read prophage_dataframe.tsv
# For each vOTU in 'Cluster' column, get list of all species in 'Classification' column for it
# Create spreadsheet with vOTUs as rows, for each column, options of what combinations of species have it
# e.g. Only E. faecium, E. faecium and E. lactis, E. faecium and E.hirae, etc.

import pandas as pd

# Read the prophage_dataframe.tsv file into a DataFrame
df = pd.read_csv('prophage_dataframe.tsv', sep='\t', low_memory=False)

# Create a new DataFrame to store vOTUs and their corresponding unique species combinations
votu_species_df = pd.DataFrame(columns=['vOTU', 'Species_Combo'])

# Create a dictionary to keep track of species for each vOTU
votu_species_dict = {}

df = df.sort_values('Classification')
# Iterate through each row in the original DataFrame
for index, row in df.iterrows():
    vOTU = row['Cluster']
    species = row['Classification']

    # Check if the vOTU already exists in the dictionary
    if vOTU in votu_species_dict:
        # Check if the species already exists for the vOTU
        if species not in votu_species_dict[vOTU]:
            # Append the species to the existing row and update the dictionary
            votu_species_df.loc[votu_species_df['vOTU'] == vOTU, 'Species_Combo'] += f',{species}'
            votu_species_dict[vOTU].append(species)
    else:
        # Create a new row for the vOTU and update the dictionary
        votu_species_dict[vOTU] = [species]
        votu_species_df = votu_species_df.append({'vOTU': vOTU, 'Species_Combo': species}, ignore_index=True)

# Write the vOTU species combinations to a new tsv file
votu_species_df.to_csv('vOTU_Species_Combo.tsv', sep='\t', index=False)

