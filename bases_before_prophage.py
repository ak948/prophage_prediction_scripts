import csv

# Input and output paths
tsv_file = "prophage_location_info.tsv"
fasta_before = "nucleotides_before.fasta"
fasta_after = "nucleotides_after.fasta"

# Open input and output files
with open(tsv_file, "r") as infile, \
     open(fasta_before, "w") as out_before, \
     open(fasta_after, "w") as out_after:

    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:
        name = row["Prophage Name"].strip()

        before_seq = row["Nucleotides Before"].strip().upper()
        after_seq = row["Nucleotides After"].strip().upper()

        # Handle 'before' sequence
        if len(before_seq) >= 1:
            out_before.write(f">{name}_before\n")
            for i in range(0, len(before_seq), 60):
                out_before.write(before_seq[i:i+60] + "\n")

        # Handle 'after' sequence
        if len(after_seq) >= 1:
            out_after.write(f">{name}_after\n")
            for i in range(0, len(after_seq), 60):
                out_after.write(after_seq[i:i+60] + "\n")

print("✅ FASTA files created:")
print(f"  • {fasta_before}")
print(f"  • {fasta_after}")
