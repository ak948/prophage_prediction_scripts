import pandas as pd

# Load the data from the TSV file
df = pd.read_csv('prophage_dataframe.tsv', sep='\t', low_memory=False)

# Extract columns related to ARG Family and Virulence Gene
arg_columns = [col for col in df.columns if col.startswith("Prophage ARG Family '")]
vir_columns = [col for col in df.columns if col.startswith("Prophage Virulence Gene '")]
def_columns = [col for col in df.columns if col.startswith("Prophage Phage Defense System Subtype '")]
ant_columns = [col for col in df.columns if col.startswith("Prophage Anti-Defense System Subtype '")]
amg_columns = [col for col in df.columns if col.startswith("Prophage AMG '")]
amg_pathway_columns = [col for col in df.columns if col.startswith("Prophage AMG Pathway '")]
amg_metabolism_columns = [col for col in df.columns if col.startswith("Prophage AMG Metabolism '")]

# Group by BiosampleID and count the number of Prophages with values greater than 0 for each gene type
biosample_summary = df.groupby('BiosampleID').apply(lambda x: (x[arg_columns] > 0).sum()).reset_index()
biosample_summary_vir = df.groupby('BiosampleID').apply(lambda x: (x[vir_columns] > 0).sum()).reset_index()
biosample_summary_def = df.groupby('BiosampleID').apply(lambda x: (x[def_columns] > 0).sum()).reset_index()
biosample_summary_ant = df.groupby('BiosampleID').apply(lambda x: (x[ant_columns] > 0).sum()).reset_index()
biosample_summary_amg = df.groupby('BiosampleID').apply(lambda x: (x[amg_columns] > 0).sum()).reset_index()
biosample_summary_amg_pathway = df.groupby('BiosampleID').apply(lambda x: (x[amg_pathway_columns] > 0).sum()).reset_index()
biosample_summary_amg_metabolism = df.groupby('BiosampleID').apply(lambda x: (x[amg_metabolism_columns] > 0).sum()).reset_index()

# Combine ARG Family and Virulence Gene summaries
biosample_summary_combined = pd.concat([biosample_summary, biosample_summary_def.drop(columns='BiosampleID'),
                                        biosample_summary_ant.drop(columns='BiosampleID'),
                                        biosample_summary_amg.drop(columns='BiosampleID'),
                                        biosample_summary_vir.drop(columns='BiosampleID'),
                                       biosample_summary_amg_pathway.drop(columns='BiosampleID'),
                                       biosample_summary_amg_metabolism.drop(columns='BiosampleID')], axis=1)

# Save the summary to a TSV file
biosample_summary_combined.to_csv('biosample_summary.tsv', sep='\t', index=False)

# Read the output TSV file containing the summary
output_file = 'biosample_summary.tsv'
output_df = pd.read_csv(output_file, sep='\t')

# Initialize a dictionary to store the counts
value_counts = {}

# Iterate over each column in the DataFrame
for column in output_df.columns[1:]:  # Exclude the 'BiosampleID' column
    # Initialize count for different values
    count_zero = 0
    count_one = 0
    count_two = 0
    count_three = 0
    count_four = 0
    count_five = 0
    count_six = 0
    count_seven = 0
    count_eight = 0
    count_nine = 0
    count_ten_or_higher = 0

    # Iterate over the values in the column
    for value in output_df[column]:
        if not value == 'No prophage genes':
            value_float = float(value)
            if int(value_float) == 0:
                count_zero += 1
            elif int(value_float) == 1:
                count_one += 1
            elif int(value_float) == 2:
                count_two += 1
            elif int(value_float) == 3:
                count_three += 1
            elif int(value_float) == 4:
                count_four += 1
            elif int(value_float) == 5:
                count_five += 1
            elif int(value_float) == 6:
                count_six += 1
            elif int(value_float) == 7:
                count_seven += 1
            elif int(value_float) == 8:
                count_eight += 1
            elif int(value_float) == 9:
                count_nine += 1
            elif int(value_float) >= 10:
                count_ten_or_higher += 1

    # Store the counts in the dictionary
    value_counts[column] = {
        '0': count_zero,
        '1': count_one,
        '2': count_two,
        '3': count_three,
        '4': count_four,
        '5': count_five,
        '6': count_six,
        '7': count_seven,
        '8': count_eight,
        '9': count_nine,
        '10+': count_ten_or_higher
    }

# Create a DataFrame from the counts
value_counts_df = pd.DataFrame(value_counts)

# Write the counts to a new TSV file
counts_output_file = 'biosample_counts_output.tsv'
value_counts_df.to_csv(counts_output_file, sep='\t', index=False)

print(f"Value counts have been saved to {counts_output_file}")







