import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import pylab
import matplotlib.pyplot as plt
from pylab import MaxNLocator
import urllib
import socket
import gzip
import time
from datetime import datetime
from Bio import Entrez


# os.chdir('Completed')
# Get list of biosampleIDs
df = pd.read_csv('genomeinfo.tsv', sep='\t')
biosampleid = df['Biosampleid']


row_num = 0
for term in biosampleid:
    search = term + "AND Enterococcus"
    retry_limit = 100  # Maximum number of retries
    retry_count = 0
    while retry_count < retry_limit:
        try:
            Entrez.email = "ak948@leicester.ac.uk"
            handle = Entrez.esearch(db="biosample", term=term)
            record = Entrez.read(handle, validate=False)
            print(record)
            # Get esummary for ID
            handle = Entrez.esummary(db="biosample", id=record["IdList"][0], rettype="docsum", retstart="xml")
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
            df.loc[df["Biosampleid"] == term, "Organism Name"] = organism_name
            df.loc[df["Biosampleid"] == term, "Collection Date"] = col_date_val
            df.loc[df["Biosampleid"] == term, "Collection Location"] = geo_loc_val
            df.loc[df["Biosampleid"] == term, "Isolation Source"] = iso_source_val
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
        except Bio.Entrez.Parser.CorruptedXMLError:
            time.sleep(10)
            retry_count += 1
            if retry_count > 95:
                time.sleep(300)
            continue
        break
df.to_csv('genomeinfo.tsv', sep='\t', index=False)