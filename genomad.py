def run_genomad(term, database, cores, genomad):
    """Run the genomad prophage prediction software against results for phageboost
	Args:
		term:genome being run through genomad, should be fasta file
	"""

    import os
    import subprocess


    directory = 'phageboost_results'  # Directory path
    base_name = f'{term}_phageboost_prophage'  # Base filename
    file_type = '.fasta'
    file_num = 1  # Starting number for renaming
    os.chdir(term)
    os.chdir(directory)

    #Run genomad on prophage genomes predicted by phageboost
    if not os.path.isfile(f'final_{term}_genomad_prophage{file_num}.fna') or os.path.isdir(f'{term}_p{file_num}'):
        for filename in os.listdir():
            if filename.endswith(file_type):
                fasta_filename = f"final_{base_name}{file_num}"
                subprocess.run(f"{genomad} end-to-end --disable-nn-classification --cleanup {fasta_filename}.fasta {fasta_filename} {database} -t {cores}",
                    shell=True) # Change database path to correct path
                #Move produced .fna results to same folder as the phageboost ones and rename them so they share a common naming scheme
                os.rename(
                    f"{fasta_filename}/{fasta_filename}_summary/{fasta_filename}_virus.fna",
                    f"final_{term}_genomad_prophage{file_num}.fna")
                file_num +=1
            else:
                None


def predicted_prophages(term):
    """Making decisions and notes on what prophages were predicted and are being kept
    Args:
        term:biosampleid for the genome the prophages were predicted in
    """

    import os
    import subprocess
    import glob

    #Deciding what prophage files to keep

    num_files = sum(1 for file in os.listdir() if file.endswith('.fasta'))
    file_num = 1

    # Iterate over the file list with a counter
    while file_num <= num_files:
        vibrant_file = f'final_{term}_vibrant_prophage{file_num}.fna'
        genomad_file = f'final_{term}_genomad_prophage{file_num}.fna'

        length_vibrant = os.stat(vibrant_file).st_size
        length_genomad = os.stat(genomad_file).st_size

        # Check if produced file is empty or not, then make sure only one prophage genome kept for each predicted
        if length_vibrant == 0 and length_genomad == 0:
            os.remove(f'{vibrant_file}')
            os.remove(f'{genomad_file}')
            with open('prophages.txt', "a") as myfile:
                myfile.write(f'\t\t\t{vibrant_file}\t{genomad_file}\tBoth empty\n')
            print(f'{vibrant_file} and {genomad_file} empty\n')
        elif length_vibrant != 0 and length_genomad == 0:
            os.remove(f'{genomad_file}')
            with open('prophages.txt', "a") as myfile:
                myfile.write(f'\t\t{vibrant_file}\t{genomad_file}\t{genomad_file} empty\n')
            print(f'{vibrant_file} kept and {genomad_file} empty\n')
        elif length_genomad != 0 and length_vibrant == 0:
            os.remove(f'{vibrant_file}')
            with open('prophages.txt', "a") as myfile:
                myfile.write(f'\t\t{genomad_file}\t{vibrant_file}\t{vibrant_file}empty\n')
            print(f'{genomad_file} kept and {vibrant_file} discarded, empty\n')
        elif length_vibrant != 0 and length_genomad !=0 and length_genomad == length_vibrant:
            os.remove(f'{genomad_file}')
            with open('prophages.txt', "a") as myfile:
                myfile.write(f'\t\t{vibrant_file}\t{genomad_file}\tFiles same size\n')
            print(f'{vibrant_file} same size\t {genomad_file} discarded\n')
        elif length_vibrant != 0 and length_genomad != 0:
            if length_vibrant > length_genomad:
                os.remove(f'{genomad_file}')
                with open('prophages.txt', "a") as myfile:
                    myfile.write(f'\t\t{vibrant_file}\t{genomad_file}\t{genomad_file} shorter\n')
                print(f'{vibrant_file} larger\t{genomad_file} discarded\n')
            if length_genomad > length_vibrant:
                os.remove(f'{vibrant_file}')
                with open('prophages.txt', "a") as myfile:
                    myfile.write(f'\t\t{genomad_file}\t{vibrant_file}\t{vibrant_file} shorter\n')
                print(f'{genomad_file} larger\t{vibrant_file} discarded\n')
        subprocess.run(f"rm -r final_{term}_phageboost_prophage{file_num} final_{term}_phageboost_prophage{file_num}.fasta",
                           shell=True)
        file_num +=1

    #Renaming files so consistent numbering scheme again

    file_extension = ".fna"  # Extension of the file you want to rename
    file_prefix = "final_" + term + "_prophage"  # Prefix of the file you want to rename
    new_file_prefix = "final_" + term + "_p"

    # Get a list of all files in the current directory
    files = os.listdir()

    # Filter and sort the files based on the prefix, term, and extension
    files = [f for f in files if f.endswith(file_extension)]
    files = sorted(files)

    # Rename the files with an increasing number
    for i, file in enumerate(files, start=1):
        new_filename = new_file_prefix + str(i) + file_extension
        os.rename(file, new_filename)
        # Make note of filename changes
        with open('prophages.txt', "a") as myfile:
            myfile.write(f'\t\t\t\t\t{file}\t{new_filename}\n')


def run_prokka(term, database, cores, prokka):
    """RUn Prokka on the results for vibrant, genomad, and phageboost to annotate the genomes
    Args:
        term:biosampleid for the genome ran through the prophage prediction softwares
    """

    import subprocess
    import os

    # Run prokka on predicted prophage genomes
    file_type = '.fna'
    base_name = f'{term}_p'
    file_num = 1
    if not os.path.isdir(f'{base_name}{file_num}') and not os.path.isdir(f'{term}_p{file_num}'):
        num_files = sum(1 for file in os.listdir() if file.endswith('.fna'))
        file_num = 1

        # Iterate over the file list with a counter, and run PROKKA on any present
        while file_num <= num_files:
            subprocess.run(f"sed '1s/.*/>{term}_p{file_num}/' final_{base_name}{file_num}{file_type} > {base_name}{file_num}{file_type}", shell=True)
            subprocess.run(
                f"{prokka} --outdir {term}_p{file_num} --prefix {term}_p{file_num} --locus {term}_p{file_num} --cpus {cores} --hmms {database} {base_name}{file_num}{file_type}",
                shell=True)
            os.chdir(f'{term}_p{file_num}')
            subprocess.run(f"rm {term}_p{file_num}.err {term}_p{file_num}.ffn {term}_p{file_num}.fsa {term}_p{file_num}.log "
                           f"{term}_p{file_num}.sqn {term}_p{file_num}.tbl {term}_p{file_num}.tsv {term}_p{file_num}.txt", shell=True)
            os.chdir('..')
            os.remove(f'{base_name}{file_num}.fna')
            file_num += 1
    else:
        os.chdir('../..')
        with open('log.txt', "a") as myfile:
            myfile.write(f'{term} prokka annotation already done\n')
        os.chdir(term)
        os.chdir('phageboost_results')

