import pandas as pd

# Load the TSV file into a DataFrame
df = pd.read_csv('results.tsv', sep='\t')

# Get the correlation matrix
corr_matrix = df.corr()

# Initialize two empty lists to store the filtered correlations
positive_correlations = []
negative_correlations = []

# Iterate over the correlation matrix and filter positive and negative correlations
for i in range(len(corr_matrix.columns)):
    for j in range(i+1, len(corr_matrix.columns)):  # Only check upper triangle of matrix
        corr_value = corr_matrix.iloc[i, j]
        if corr_value > 0.5:  # Positive correlations
            positive_correlations.append([corr_matrix.columns[i], corr_matrix.columns[j], corr_value])
        elif corr_value < -0.5:  # Negative correlations
            negative_correlations.append([corr_matrix.columns[i], corr_matrix.columns[j], corr_value])

# Convert the lists to DataFrames
positive_df = pd.DataFrame(positive_correlations, columns=['First Column', 'Second Column', 'Correlation Value'])
negative_df = pd.DataFrame(negative_correlations, columns=['First Column', 'Second Column', 'Correlation Value'])

positive_df.to_csv(f'results_corr_positive.tsv', sep='\t', index=False)
negative_df.to_csv(f'results_corr_negative.tsv', sep='\t', index=False)
