import pandas as pd

# Load the TSV file
df = pd.read_csv('../../bacteria_amr_genes.tsv', sep='\t')
# Create a dictionary where each drug class (antibiotic) maps to a set of genes
drug_to_genes = {}
# Iterate through each row in the dataframe
for _, row in df.iterrows():
    gene_family = row['AMR Gene Family']
    if pd.isna(row['Drug Class']):
        continue
    drugs = [drug.strip() for drug in row['Drug Class'].split(';')]
    for drug in drugs:
        drug_to_genes.setdefault(drug, set()).add(gene_family)
# Build the final output DataFrame
all_drugs = sorted(drug_to_genes.keys())
max_rows = max(len(genes) for genes in drug_to_genes.values())
output_data = {drug: list(genes) + [''] * (max_rows - len(genes)) for drug, genes in drug_to_genes.items()}
output_df = pd.DataFrame(output_data)
# Save to TSV
output_df.to_csv('genes_by_drug_class.tsv', sep='\t', index=False)

# Load the TSV file
df = pd.read_csv('../../bacteria_phage_defenses.tsv', sep='\t')
# Create a dictionary where each defense type maps to a set of genes (defense subtypes)
type_to_genes = {}
# Iterate through each row in the dataframe
for _, row in df.iterrows():
    defense_subtype = row['Defense Subtype']
    if pd.isna(row['Defense Type']):
        continue
    defense = row['Defense Type']
    type_to_genes.setdefault(defense, set()).add(defense_subtype)
# Build the final output DataFrame
all_defenses = sorted(type_to_genes.keys())
max_rows = max(len(genes) for genes in type_to_genes.values())
output_data = {defense: list(genes) + [''] * (max_rows - len(genes)) for defense, genes in type_to_genes.items()}
output_df = pd.DataFrame(output_data)
# Save to TSV
output_df.to_csv('defenses_by_defense_type.tsv', sep='\t', index=False)

