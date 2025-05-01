import pandas as pd
import re

infile = "new.tsv"
outfile = "cleaned_my_clusters.tsv"

df = pd.read_csv('removed.txt', sep='\t')
delete_list = df['Biosampleid']
with open(infile) as fin, open(outfile, "w+") as fout:
    for line in fin:
        for word in delete_list:
            line = re.sub(r'\b{}\b[,\"\t]'.format(re.escape(word)), "", line)
        fout.write(line)

