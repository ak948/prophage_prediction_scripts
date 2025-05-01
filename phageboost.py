def get_url(term, access, email, download=True):
    """Download fasta file for a given search term.
	Args:
		term: search term, usually biosampleid
		access: Accession number
		download: whether to download the results
		path:folder to save to
	"""

    import os
    import urllib
    import socket
    import gzip
    import shutil
    import subprocess
    import pandas as pd
    import time
    from datetime import datetime
    from Bio import Entrez

    # Set email for Entrez
    Entrez.email = email
    # Check if folder already made
    if not os.path.isdir(term) and not os.path.isfile(f'{term}.tar.gz'):
        df = pd.read_csv('genomeinfo.tsv', delimiter='\t')
        filtered_df = df[~df['Biosampleid'].str.contains(term)]
        filtered_df.to_csv('genomeinfo.tsv', sep='\t', index=False)
        retry_limit = 100  # Maximum number of retries
        retry_count = 0
        while retry_count < retry_limit:
            try:
                # Get ID for BioSampleID
                handle = Entrez.esearch(db="assembly", term=access)
                record = Entrez.read(handle, validate=False)
                # Get esummary for ID
                handle = Entrez.efetch(db="assembly", id=record["IdList"], rettype="docsum", retmode="xml")
                record = Entrez.read(handle, validate=False)
                # Get last update date
                date = record['DocumentSummarySet']['DocumentSummary'][0]['LastUpdateDate']
                os.mkdir(term)
                with open(f'{term}/info.txt', "a") as myfile:
                    myfile.write(f'{date}\t{access}\n')
                # get ftp link
                url = record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                label = os.path.basename(url)
                # get the fasta link - change to get other formats
                link = os.path.join(url, label + '_genomic.fna.gz')
                print(link)
                if download == True:
                    # download link
                    timeout = 10
                    socket.setdefaulttimeout(timeout)
                    while True:
                        try:
                            urllib.request.urlretrieve(link, f'{label}.fna.gz')
                        except urllib.error.URLError as e:
                            time.sleep(1)
                            continue
                        break
                    # Rename fasta file to match BioSampleID and put in folder
                    os.rename(f'{label}.fna.gz', f'{term}/{term}.fna.gz')
                    import gzip
                    import shutil
                    with gzip.open(f'{term}/{term}.fna.gz', 'rb') as f_in:
                        with open(f'{term}/{term}.fna', 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(f'{term}/{term}.fna.gz')

                else:
                    os.chdir('..')
                    with open('log.txt', "a") as myfile:
                        myfile.write(f'{term} genome already downloaded\n')
                    os.chdir(term)
            except urllib.error.URLError as e:
                time.sleep(10)
                retry_count += 1
                if retry_count > 95:
                    time.sleep(300)
                continue
            except IndexError as e:
                time.sleep(10)
                retry_count += 1
                if retry_count > 95:
                    time.sleep(300)
                continue
            except RuntimeError as e:
                time.sleep(10)
                retry_count += 1
                if retry_count > 95:
                    time.sleep(300)
                continue
            except socket.timeout:
                time.sleep(10)
                retry_count += 1
                if retry_count > 95:
                    time.sleep(300)
                continue
            break


    else:
        if os.path.isfile(f'{term}.tar.gz'):
            subprocess.run(f"tar -zxf {term}.tar.gz {term}", shell=True)
            os.remove(f'{term}.tar.gz')

        retry_limit = 100  # Maximum number of retries
        retry_count = 0
        while retry_count < retry_limit:
            try:
                # Read date from file
                with open(f'{term}/info.txt', 'r') as file:
                    line = file.readline().strip()  # Read the first line
                    date1 = line.split('\t')[0]  # Extract the date part
                    accession1 = line.split('\t')[1]  # Extract accession number

                handle = Entrez.esearch(db="assembly", term=access)
                record = Entrez.read(handle, validate=False)
                # Get esummary for ID
                handle = Entrez.efetch(db="assembly", id=record["IdList"], rettype="docsum", retmode="xml")
                record = Entrez.read(handle, validate=False)
                date = record['DocumentSummarySet']['DocumentSummary'][0]['LastUpdateDate']

                # Convert both dates to datetime objects
                date1_obj = datetime.strptime(date1, "%Y/%m/%d %H:%M")
                date2_obj = datetime.strptime(date, "%Y/%m/%d %H:%M")

                if not access == accession1:
                    # Compare the dates
                    if date1_obj > date2_obj:  # If new biosampleID older, no changes
                        with open('repeats.txt', 'a') as myfile:
                            myfile.write(f'{term}\t{date}\t{access}\nOlder than {accession1}, no changes')
                        print(f"Repeat biosampleID, new accession, {access} older than {accession1}, no changes")
                        subprocess.run(f"tar -zcf {term}.tar.gz {term}", shell=True)
                        shutil.rmtree(term)
                    elif date1_obj == date2_obj:
                        with open('repeats.txt', 'a') as myfile:
                            myfile.write(f'{term}\t{date}\t{access}\nSame age as {accession1}, no changes')
                        print(f"Repeat biosampleID, new accession, {access} same age as {accession1}, no changes")
                        subprocess.run(f"tar -zcf {term}.tar.gz {term}", shell=True)
                        shutil.rmtree(term)
                    elif date2_obj > date1_obj:  # If new biosampleID newer, make note and remove original directory to repeat
                        with open('repeats.txt', 'a') as myfile:
                            myfile.write(f'{term}\t{date1}\t{accession1}\nOlder than {access}, redoing')
                        print(
                            f"Repeat biosampleID, old accession {accession1} older than {access}, redoing for new accession number")
                        subprocess.run(f"rm -r {term}", shell=True)
                        # Remove old genome info
                        df = pd.read_csv('genomeinfo.tsv', delimiter='\t')
                        filtered_df = df[~df['Biosampleid'].str.contains(term)]
                        filtered_df.to_csv('genomeinfo.tsv', sep='\t', index=False)
                        Entrez.email = email
                        # Get ID for BioSampleID
                        handle = Entrez.esearch(db="assembly", term=access)
                        record = Entrez.read(handle, validate=False)
                        # Get esummary for ID
                        handle = Entrez.efetch(db="assembly", id=record["IdList"], rettype="docsum", retmode="xml")
                        record = Entrez.read(handle, validate=False)
                        # Get update date and accession number
                        date = record['DocumentSummarySet']['DocumentSummary'][0]['LastUpdateDate']
                        os.mkdir(term)
                        with open(f'{term}/info.txt', "a") as myfile:
                            myfile.write(f'{date}\t{access}\n')
                        # get ftp link
                        url = record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
                        label = os.path.basename(url)
                        # get the fasta link - change to get other formats
                        link = os.path.join(url, label + '_genomic.fna.gz')
                        print(link)
                        if download == True:
                            # download link
                            urllib.request.urlretrieve(link, f'{label}.fna.gz')
                            # Rename fasta file to match BioSampleID and put in folder
                            os.rename(f'{label}.fna.gz', f'{term}/{term}.fna.gz')
                            import gzip
                            import shutil
                            with gzip.open(f'{term}/{term}.fna.gz', 'rb') as f_in:
                                with open(f'{term}/{term}.fna', 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                            os.remove(f'{term}/{term}.fna.gz')
            except urllib.error.URLError as e:
                time.sleep(10)
                retry_count += 1
                if retry_count > 95:
                    time.sleep(300)
                continue
            except IndexError as e:
                time.sleep(10)
                retry_count += 1
                if retry_count > 95:
                    time.sleep(300)
                continue
            except RuntimeError as e:
                time.sleep(10)
                retry_count += 1
                if retry_count > 95:
                    time.sleep(300)
                continue
            except socket.timeout:
                time.sleep(10)
                retry_count += 1
                if retry_count > 95:
                    time.sleep(300)
                continue
            break


def get_assembly_info(term, email):
    """Record additional info for a given search term.
    Args:
		term:search term, usually biosampleid
	"""
    import csv
    import os
    from Bio import Entrez
    from datetime import datetime
    import socket
    import time
    import urllib

    if not os.path.isfile(f'{term}.tar.gz'):
        os.chdir(term)
        if not os.path.isdir('phageboost_results'):
            retry_limit = 100  # Maximum number of retries
            retry_count = 0
            while retry_count < retry_limit:
                try:
                    # Set email for Entrez
                    Entrez.email = email
                    # Get ID for BioSampleID
                    handle = Entrez.esearch(db="assembly", term=term)
                    record = Entrez.read(handle, validate=False)
                    biosampleid = record["IdList"]
                    # Get esummary for ID
                    handle = Entrez.efetch(db="assembly", id=biosampleid, rettype="docsum", retmode="xml")
                    record = Entrez.read(handle, validate=False)
                    # get organism, submission date, last updated date and submitter organisation
                    organism = record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
                    submissiondate = record['DocumentSummarySet']['DocumentSummary'][0]['SubmissionDate']
                    lastupdatedate = record['DocumentSummarySet']['DocumentSummary'][0]['LastUpdateDate']
                    submitter = record['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']
                    # Get ID for BioSampleID
                    try:
                        handle = Entrez.esearch(db="biosample", term=term)
                        record = Entrez.read(handle, validate=False)
                        # Get esummary for ID
                        handle = Entrez.esummary(db="biosample", id=record["IdList"][0], rettype="docsum",
                                                 retstart="xml")
                        record = Entrez.read(handle, validate=False)
                        # Parse the XML string
                        attributes = record['DocumentSummarySet']['DocumentSummary'][0]['SampleData']
                        import xml.etree.ElementTree as ET
                        root = ET.fromstring(attributes)

                        # Extract organism name, collection date, collecion location, isolation source
                        OrganismName = root.find('.//OrganismName')
                        organism_name = OrganismName.text if OrganismName is not None else "Missing"
                        collection_date = root.find(
                            ".//*[@attribute_name='collection_date']")  # Find the attribute with attribute_name = "collection_date"
                        col_date_val = collection_date.text if collection_date is not None else "Missing"  # Extract the attribute value or set to "Missing" if not found
                        geo_loc = root.find(".//*[@attribute_name='geo_loc_name']")
                        geo_loc_val = geo_loc.text if geo_loc is not None else "Missing"
                        isolation_source = root.find(".//*[@attribute_name='isolation_source']")
                        iso_source_val = isolation_source.text if isolation_source is not None else "Missing"

                        # Add information to .tsv file
                        os.chdir('..')
                        with open('genomeinfo.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow(
                                [f'{term}', f'{organism}', f'{submissiondate}', f'{lastupdatedate}', f'{submitter}',
                                 f'{organism_name}',
                                 f'{col_date_val}', f'{geo_loc_val}', f'{iso_source_val}'])
                        os.chdir(term)
                    except:
                        # Add information to .tsv file
                        os.chdir('..')
                        with open('genomeinfo.tsv', 'a') as out_file:
                            tsv_writer = csv.writer(out_file, delimiter='\t')
                            tsv_writer.writerow(
                                [f'{term}', f'{organism}', f'{submissiondate}', f'{lastupdatedate}', f'{submitter}',
                                 'Missing',
                                 'Missing', 'Missing', 'Missing'])
                        os.chdir(term)
                except urllib.error.URLError as e:
                    time.sleep(10)
                    retry_count += 1
                    if retry_count > 95:
                        time.sleep(300)
                    continue
                except IndexError as e:
                    time.sleep(10)
                    retry_count += 1
                    if retry_count > 95:
                        time.sleep(300)
                    continue
                except RuntimeError as e:
                    time.sleep(10)
                    retry_count += 1
                    if retry_count > 95:
                        time.sleep(300)
                    continue
                except socket.timeout:
                    time.sleep(10)
                    retry_count += 1
                    if retry_count > 95:
                        time.sleep(300)
                    continue
                break


def run_phageboost(term, cores, phageboost):
    """Run the phageboost prophage prediction software
	Args:
		term:genome being run through phageboost, should be fasta file
	"""

    import subprocess
    import os

    if not os.path.isfile(f'{term}.tar.gz'):
        # Run phageboost on any .fna files in the current folder, then print phageboost finished once done
        if not os.path.isdir('phageboost_results'):
            os.mkdir('phageboost_results')
            subprocess.run(f"{phageboost} -f *.fna -o phageboost_results -j {cores}", shell=True)
            print("PhageBoost finished")

            # Renaming phage fasta files produced so they share a common naming scheme and are numbered
            directory = 'phageboost_results'  # Directory path
            base_name = f'{term}_phageboost_prophage'  # New base name
            file_type = '.fasta'  # File type to change the names of

            file_num = 1  # Starting number for renaming

            with open(f'{directory}/prophages.txt', "a") as myfile:
                myfile.write(
                    f'Original Filename\t1st Change\tKeep\tDiscard 1\tDiscard 2\tComment\t2nd Name\tFinal Name\n')

            for filename in os.listdir(directory):
                if filename.endswith(file_type):
                    new_filename = f"final_{base_name}{file_num}{file_type}"
                    file_num += 1
                    os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
                    # Make file to keep note of prophage name changes and what are kept and deleted
                    with open(f'{directory}/prophages.txt', "a") as myfile:
                        myfile.write(f'{filename}\t{new_filename}\n')
            subprocess.run("rm -r *.fna", shell=True)


def also_vibrant(term, cores, vibrant):
    """Runt the VIBRANT prophage prediction software against the Phageboost results
    Args:
        term:biosampleID of the predicted prophages being ran through VIBRANT
    """

    import subprocess
    import os

    if not os.path.isfile(f'{term}.tar.gz'):
        # Information for fasta files
        directory = 'phageboost_results'  # Directory path
        base_name = f'{term}_phageboost_prophage'  # Base filename
        file_type = '.fasta'
        file_num = 1  # Starting number for renaming
        os.chdir(directory)

        # Running vibrant on each fasta file produced by phageboost
        if not os.path.isfile(f'final_{term}_vibrant_prophage{file_num}.fna') or os.path.isdir(
                f'{term}_prophage{file_num}'):
            for filename in os.listdir():
                if filename.endswith(file_type):
                    fasta_filename = f"final_{base_name}{file_num}"
                    subprocess.run(
                        f"python3 {vibrant} -i {fasta_filename}.fasta -folder {fasta_filename} -t {cores}",
                        shell=True)  # Change path to VIBRANT to correct path
                    # Move produced .fna file to same folder as the phageboost results and rename them so there is a common naming scheme
                    os.rename(
                        f"{fasta_filename}/VIBRANT_{fasta_filename}/VIBRANT_phages_{fasta_filename}/{fasta_filename}.phages_combined.fna",
                        f"final_{term}_vibrant_prophage{file_num}.fna")
                    file_num += 1
                else:
                    None
