import os
import csv

with open('removed.txt', 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    # Iterate over each row in the file
    for row in reader:
        id = row
        for filename in os.listdir():
            if filename.startswith(id):
                subprocess.run(f'mv {row}.fna filtered', shell=True)