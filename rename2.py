import os
import csv
import pandas as pd
import subprocess
import shutil

os.chdir('Completed')
#Get list of biosampleIDs
df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid'].loc[0:600]

# Iterate over the ids in biosampleid
for id in biosampleid:
    file_num = 1
    if os.path.isfile(f'{id}.tar.gz'):
        subprocess.run(f"tar -zxf {id}.tar.gz {id}", shell=True)
        os.remove(f'{id}.tar.gz')
        os.chdir(f'{id}/phageboost_results')

    # Iterate over possible folder names
    while True:
        folder_name = f"{id}_prophage{file_num}"
        if os.path.exists(folder_name):
            path = f'/home/ak948/prophage_prediction/Completed/{id}/phageboost_results'
            old_path = os.path.join(path, folder_name)
            new_folder_name = f"{id}_p{file_num}"
            new_path = os.path.join(path, new_folder_name)
            os.rename(old_path, new_path)

        folder_name = f"{id}_p{file_num}"

        # Change directory to go into the folder
        if os.path.exists(folder_name):
            os.chdir(folder_name)
            path = f'/home/ak948/prophage_prediction/Completed/{id}/phageboost_results/{folder_name}'
            # Iterate over the files in the folder
            for file in os.listdir():
                old_name = f"{id}_prophage{file_num}.fna"
                if os.path.isfile(old_name):
                    old_path = os.path.join(path, old_name)
                    new_name = f"{id}_p{file_num}.fna"
                    new_path = os.path.join(path, new_name)
                    os.rename(old_path, new_path)
                old_name = f"{id}_prophage{file_num}.faa"
                if os.path.isfile(old_name):
                    old_path = os.path.join(path, old_name)
                    new_name = f"{id}_p{file_num}.faa"
                    new_path = os.path.join(path, new_name)
                    os.rename(old_path, new_path)
                old_name = f"{id}_prophage{file_num}.gbk"
                if os.path.isfile(old_name):
                    old_path = os.path.join(path, old_name)
                    new_name = f"{id}_p{file_num}.gbk"
                    new_path = os.path.join(path, new_name)
                    os.rename(old_path, new_path)
                old_name = f"{id}_prophage{file_num}.gff"
                if os.path.isfile(old_name):
                    old_path = os.path.join(path, old_name)
                    new_name = f"{id}_p{file_num}.gff"
                    new_path = os.path.join(path, new_name)
                    os.rename(old_path, new_path)

            # Change directory back to previous directory
            os.chdir("..")

        # Increment file_num for the next iteration
        file_num += 1

        # Break the loop if no folder exists for the next file_num
        if not os.path.exists(f"{id}_prophage{file_num}") and not os.path.exists(f"{id}_p{file_num}"):
            break
    os.chdir('/home/ak948/prophage_prediction/Completed')
    subprocess.run(f"tar -zcvf {id}.tar.gz {id}", shell=True)
    shutil.rmtree(id)
