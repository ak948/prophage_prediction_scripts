import csv
from collections import defaultdict, Counter
from itertools import combinations
from scipy.stats import fisher_exact
import numpy as np

input_file = "prophage_dataframe.tsv"
output_file = "votu_species_cooccurrence.tsv"

# Step 1: Map BiosampleID to set of vOTUs
biosample_to_votus = defaultdict(set)

with open(input_file, newline='') as infile:
    reader = csv.DictReader(infile, delimiter="\t")
    for row in reader:
        biosample = row["BiosampleID"].strip()
        votu = row["Enterococcus Prophage Species"].strip()
        if biosample and votu:
            biosample_to_votus[biosample].add(votu)

# Count total occurrences of each vOTU
votu_counts = Counter()
for votus in biosample_to_votus.values():
    for votu in votus:
        votu_counts[votu] += 1

# Filter to only include vOTUs with at least 20 prophages
valid_votus = {votu for votu, count in votu_counts.items() if count >= 20}

# Step 2: Count co-occurrence between vOTU pairs (only using valid vOTUs)
cooccurrence_counts = Counter()
for votus in biosample_to_votus.values():
    filtered_votus = [v for v in votus if v in valid_votus]
    for v1, v2 in combinations(sorted(filtered_votus), 2):
        cooccurrence_counts[(v1, v2)] += 1

# Step 3: Prepare to collect results and p-values for BH correction
results = []
p_values = []

for (v1, v2), together in cooccurrence_counts.items():
    total_v1 = votu_counts[v1]
    total_v2 = votu_counts[v2]

    not_with_v2 = total_v1 - together
    not_with_v1 = total_v2 - together

    # Build contingency table excluding 'neither' samples:
    #        vOTU2 present  vOTU2 absent
    # vOTU1 present    together     not_with_v2
    # vOTU1 absent     not_with_v1  0 (excluded)
    contingency_table = [
        [together, not_with_v2],
        [not_with_v1, 0]  # zero because isolates with neither are excluded
    ]

    # Fisher's exact test (two-sided)
    try:
        oddsratio, p = fisher_exact(contingency_table, alternative='two-sided')
    except Exception:
        p = 1.0
        oddsratio = float('nan')

    jaccardindex = together / (together + not_with_v1 + not_with_v2)
    togetherness_ratio = (together / (not_with_v1 + together)) / (together / (not_with_v2 + together))

    p_values.append(p)

    results.append({
        "v1": v1,
        "v2": v2,
        "together": together,
        "not_with_v2": not_with_v2,
        "not_with_v1": not_with_v1,
        "pct_v1": (together / total_v1) * 100 if total_v1 else 0,
        "pct_v2": (together / total_v2) * 100 if total_v2 else 0,
        "p": p,
        "odds_ratio": oddsratio,
        "jaccard_index": jaccardindex,
        "togethernessratio": togetherness_ratio
    })

# Step 4: Perform Benjamini-Hochberg correction on p-values
def benjamini_hochberg(pvals, alpha=0.05):
    pvals = np.array(pvals)
    n = len(pvals)
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]
    adjusted = np.empty(n, dtype=float)
    cummin = 1.0
    for i in reversed(range(n)):
        rank = i + 1
        bh_value = (n / rank) * sorted_pvals[i]
        cummin = min(cummin, bh_value)
        adjusted[i] = cummin
    adjusted_pvals = np.empty(n, dtype=float)
    adjusted_pvals[sorted_indices] = adjusted
    adjusted_pvals = np.minimum(adjusted_pvals, 1.0)
    return adjusted_pvals

adjusted_pvals = benjamini_hochberg(p_values)

# Step 5: Write results with BH adjusted p-value and significance call
with open(output_file, "w", newline='') as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    writer.writerow([
        "vOTU 1", "vOTU 2", "Times Found Together",
        "vOTU 1 Not With vOTU 2", "vOTU 2 Not With vOTU 1",
        "% Times Found Together (of vOTU 1)", "% Times Found Together (of vOTU 2)",
        "Fisher's Exact p-value", "BH-adjusted p-value", "Significantly Found Together (BH)",
        "Odds Ratio", "Odds Ratio Association", "Jaccard Index", "Ratio of Found Together for vOTU 1 vs vOTU 2"
    ])

    for res, adj_p in zip(results, adjusted_pvals):
        significant = "Yes" if adj_p <= 0.05 else "No"

        if res["odds_ratio"] >= 1.25:
            odds = "Positive Association"
        elif res["odds_ratio"] <= 0.75:
            odds = "Negative Association"
        else:
            odds = "No Association"

        writer.writerow([
            res["v1"], res["v2"], res["together"],
            res["not_with_v2"], res["not_with_v1"],
            f"{res['pct_v1']:.2f}", f"{res['pct_v2']:.2f}",
            f"{res['p']:.4e}", f"{adj_p:.4e}", significant,
            f"{res['odds_ratio']:.3f}", odds, res["jaccard_index"], res["togethernessratio"]
        ])

isolates_with_one_votu = sum(
    1 for vset in biosample_to_votus.values() if len([v for v in vset if v in valid_votus]) == 1
)

print(f"✅ Saved co-occurrence table to {output_file}")
print(f"🔹 Isolates with only one vOTU (≥20 prophages): {isolates_with_one_votu}")
