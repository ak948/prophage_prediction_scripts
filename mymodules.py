def get_url(term, download=True):
    """Download fasta file for a given search term.
	Args:
		term: search term, usually biosampleid
		download: whether to download the results
		path:folder to save to
	"""

    import os
    import urllib
    from Bio import Entrez

    # Set email for Entrez
    Entrez.email = "ak948@leicester.ac.uk"
    # Check if folder already made
    if not os.path.isdir(term):
        # Get ID for BioSampleID
        handle = Entrez.esearch(db="assembly", term=term)
        record = Entrez.read(handle)
        # Get esummary for ID
        handle = Entrez.efetch(db="assembly", id=record["IdList"], rettype="docsum", retmode="xml")
        record = Entrez.read(handle)
        # get ftp link
        url = record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        label = os.path.basename(url)
        # get the fasta link - change to get other formats
        link = os.path.join(url, label + '_genomic.fna.gz')
        print(link)
        if download == True:
            # make folder
            os.mkdir(term)
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
            with open(f'{term}/{term}_log.txt', "w") as myfile:
                myfile.write(f'Genome downloaded successfully\n')
        else:
            os.chdir('..')
            with open('log.txt', "a") as myfile:
                myfile.write(f'{term} genome already downloaded\n')
            os.chdir(term)

def get_assembly_info(term):
    """Record additional info for a given search term.
    Args:
		term:search term, usually biosampleid
	"""
    import csv
    import os
    if os.path.isdir(term):
        os.chdir(term)
        if not os.path.isfile(f'{term}.txt'):
            from Bio import Entrez
            # Set email for Entrez
            Entrez.email = "ak948@leicester.ac.uk"
            # Get ID for BioSampleID
            handle = Entrez.esearch(db="assembly", term=term)
            record = Entrez.read(handle)
            biosampleid = record["IdList"]
            # Get esummary for ID
            handle = Entrez.efetch(db="assembly", id=biosampleid, rettype="docsum", retmode="xml")
            record = Entrez.read(handle)
            # get organism, submission date, last updated date and submitter organisation
            organism = record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
            submissiondate = record['DocumentSummarySet']['DocumentSummary'][0]['SubmissionDate']
            lastupdatedate = record['DocumentSummarySet']['DocumentSummary'][0]['LastUpdateDate']
            submitter = record['DocumentSummarySet']['DocumentSummary'][0]['SubmitterOrganization']
            #Get ID for BioSampleID
            handle = Entrez.esearch(db="biosample", term=term)
            record = Entrez.read(handle)
            # Get esummary for ID
            handle = Entrez.esummary(db="biosample", id=record["IdList"], rettype="docsum", retstart="xml")
            record = Entrez.read(handle)
            # Parse the XML string
            attributes = record['DocumentSummarySet']['DocumentSummary'][0]['SampleData']
            import xml.etree.ElementTree as ET
            root = ET.fromstring(attributes)

            # Extract organism name, collection date, collecion location, isolation source
            organism_name = root.find('.//OrganismName').text
            collection_date = root.find(
                ".//*[@attribute_name='collection_date']")  # Find the attribute with attribute_name = "collection_date"
            col_date_val = collection_date.text if collection_date is not None else "Missing"  # Extract the attribute value or set to "Missing" if not found
            geo_loc = root.find(".//*[@attribute_name='geo_loc_name']")
            geo_loc_val = geo_loc.text if geo_loc is not None else "Missing"
            isolation_source = root.find(".//*[@attribute_name='isolation_source']")
            iso_source_val = isolation_source.text if isolation_source is not None else "Missing"

            #Add information to .tsv file
            os.chdir('..')
            with open('genomeinfo.tsv', 'a') as out_file:
                tsv_writer = csv.writer(out_file, delimiter='\t')
                tsv_writer.writerow([f'{biosampleid}', f'{organism}', f'{submissiondate}', f'{lastupdatedate}', f'{submitter}', f'{organism_name}',
                     f'{col_date_val}', f'{geo_loc_val}', f'{iso_source_val}'])
            os.chdir(term)

        else:
            os.chdir('..')
            with open('log.txt', "a") as myfile:
                myfile.write(f'{term} info already added\n')
            os.chdir(term)


def run_phispy(term):
    """Run the PhiSpy prophage prediction software
	Args:
			term:genome being run through PhiSpy, should be fasta file
	"""

    import subprocess
    import os
    import glob

    #Run Prokka on .fna file to get annotated genbank file
    if not glob.glob('PROKKA_*'):
        subprocess.run("prokka --proteins FASTA *.fna", shell=True)
    else:
        print("PROKKA already done")

    #Run PhiSpy on .gbk file produced by Prokka
    if not os.path.isdir('phispy_results'):
        subprocess.run("PhiSpy.py PROKKA*/*.gbk -o phispy_results --output_choice 7", shell=True)
        print("PhiSpy finished")
    else:
        os.chdir('..')
        with open('log.txt', "a") as myfile:
            myfile.write(f'{term} phispy prediction already done\n')
        os.chdir(term)

def run_genomad(term):
    """Run the genomad prophage prediction software
	Args:
		term:genome being run through genomad, should be fasta file
	"""

    import subprocess
    import os

    #Run genomad on any .fna files in the current folder, then print genomad finished once done
    if not os.path.isdir('genomad_results'):
        subprocess.run("genomad end-to-end --disable-nn-classification --cleanup *.fna genomad_results /data/DB/genomad_db", shell=True)
        print("genomad finished")
    else:
        os.chdir('..')
        with open('log.txt', "a") as myfile:
            myfile.write(f'{term} genomad prediction already done\n')
        os.chdir(term)


def run_vibrant(term):
    """Run the VIBRANT prophage prediction software
	Args:
		term:genome being run through VIBRANT, should be fasta file
	"""

    import subprocess
    import os

    # Run VIBRANT on any .fna files in the current folder, then print VIBRANT finished once done
    if not os.path.isdir('vibrant_results'):
        subprocess.run("python3 /data/DB/VIBRANT/VIBRANT_run.py -i *.fna -folder vibrant_results -t 16", shell=True)
        print("VIBRANT finished")

    else:
        os.chdir('..')
        with open('log.txt', "a") as myfile:
            myfile.write(f'{term} vibrant prediction already done\n')
        os.chdir(term)