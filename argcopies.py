import pandas as pd

# Read the input TSV file
input_file = 'results.tsv'
df = pd.read_csv(input_file, sep='\t')

# Initialize a list to hold results for the new file
results = []
column_headers = ['BiosampleID']

# Set to track which columns we need in the header (to avoid duplication)
included_columns = set()

# Iterate over each row in the DataFrame
for idx, row in df.iterrows():
    biosample_id = row['BiosampleID']
    differences = []

    # Flag to track if we have added any valid differences for this row
    has_prophage_gene = False

    # Iterate over columns to find ARG Family and Virulence Gene columns
    for column in df.columns:
        if column.startswith("Bacteria ARG Family '") or column.startswith("Bacteria Virulence Gene '") or column.startswith("Bacteria Phage Defense System Subtype '"):
            # Extract the ARG Family or Virulence Gene name from the column header
            gene_name = column.split("'")[1]  # Extract name within quotes

            # Determine if the column is for an ARG Family or Virulence Gene
            if column.startswith('Bacteria ARG Family'):
                gene_type = "ARG Family"
            elif column.startswith('Bacteria Virulence Gene'):
                gene_type = "Virulence Gene"
            else:
                gene_type = "Phage Defense System Subtype"

            # Find the corresponding Prophage column
            prophage_column = f"Prophage {gene_type} '{gene_name}' Count"

            # Check if the Prophage column exists in the DataFrame
            if prophage_column in df.columns:
                arg_value = row[column]
                prophage_value = row[prophage_column]

                # If the Prophage value is 0, skip the difference for this gene
                if prophage_value > 0:
                    difference = arg_value - prophage_value
                    differences.append(difference)

                    # Add the Prophage column to the result headers if it's the first time we encounter it
                    if f"{gene_type}: {gene_name} Diff" not in included_columns:
                        column_headers.append(f"{gene_type}: {gene_name} Diff")
                        included_columns.add(f"{gene_type}: {gene_name} Diff")

                    # Mark that we've found a Prophage gene with value > 0
                    has_prophage_gene = True
                else:
                    difference = "No prophage genes"
                    differences.append(difference)

                    # Add the Prophage column to the result headers if it's the first time we encounter it
                    if f"{gene_type}: {gene_name} Diff" not in included_columns:
                        column_headers.append(f"{gene_type}: {gene_name} Diff")
                        included_columns.add(f"{gene_type}: {gene_name} Diff")

    # Only add the row if we have valid data (i.e., some Prophage values were > 0)
    if any(value != "No prophage genes" for value in differences):
        # Ensure the number of columns matches the number of headers by padding the differences list
        results.append([biosample_id] + differences)

# Create the output DataFrame
results_df = pd.DataFrame(results, columns=column_headers)

# Write the results to a new TSV file
output_file = 'differences_output.tsv'
results_df.to_csv(output_file, sep='\t', index=False)

print(f"Differences have been saved to {output_file}")

# Read the output TSV file containing the differences
output_file = 'differences_output.tsv'
output_df = pd.read_csv(output_file, sep='\t')

# Initialize a dictionary to store the counts
value_counts = {}

# Iterate over each column in the DataFrame
for column in output_df.columns[1:]:  # Exclude the 'BiosampleID' column
    # Initialize count for different values
    count_zero_or_lower = 0
    count_one = 0
    count_two = 0
    count_three = 0
    count_four = 0
    count_five = 0
    count_six_or_greater = 0
    print(column)
    # Iterate over the values in the column
    for value in output_df[column]:
        if not value == 'No prophage genes':
            value_float = float(value)
            if int(value_float) < 0:
                count_zero_or_lower += 1
            elif int(value_float) == 0:
                count_zero_or_lower += 1
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
            elif int(value_float) >= 6:
                count_six_or_greater += 1

    # Store the counts in the dictionary
    value_counts[column] = {
        '0': count_zero_or_lower,
        '1': count_one,
        '2': count_two,
        '3': count_three,
        '4': count_four,
        '5': count_five,
        '6 or greater': count_six_or_greater
    }

# Create a DataFrame from the counts
value_counts_df = pd.DataFrame(value_counts)

# Write the counts to a new TSV file
counts_output_file = 'value_counts_output.tsv'
value_counts_df.to_csv(counts_output_file, sep='\t', index=False)

print(f"Value counts have been saved to {counts_output_file}")



