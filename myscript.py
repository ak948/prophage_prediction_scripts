import pandas as pd
import os
import subprocess
import csv

# read excel file, setting sheet name and first row as header
df = pd.read_excel('Enterococcus.xls', sheet_name='prokaryotes', header=0)
biosampleid = df['BioSample'].head(5)
with open('genomeinfo.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['Biosampleid', 'Organism', 'Submission Date', 'Last Update Date', 'Submitter', 'Organism Name', 'Collection Date', 'Collection Location', 'Isolation Source'])
with open('log.txt', 'a') as myfile:
    myfile.write('Repeats')

from mymodules import *

# download fasta genome for biosampleid, get info for biosampleid and add to file
for id in biosampleid:
    links = get_url(id, download=True)
    assemblyinfo = get_assembly_info(id)
    runphispy = run_phispy(id)
    rungenomad = run_genomad(id)
    #cleanup
    subprocess.run("rm -r *.fna PROK*", shell=True)
    os.chdir('..')

#Once finished, create file saying so
with open(f'myscript.txt', "w") as myfile:
                myfile.write(f'Script finished')
