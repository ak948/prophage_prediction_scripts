# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import csv

data = pd.read_csv("Isolation_source.tsv", sep="\t")
animal_source = data['Animal'].dropna()
environmental_source = data['Environmental Sample'].dropna()
human_source = data['Human'].dropna()
plant_source = data['Plants and Products'].dropna()
other_source = data['Other'].dropna()
unknown_source = data['Not Provided'].dropna()
#Set the box plot data
columns = [animal_source, environmental_source, human_source, plant_source, other_source, unknown_source]
#Generate the figure
fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)
# Choose data used, add notches and colours
ax.boxplot(columns, notch=True, patch_artist=True)
#Provide x-axis names
plt.xticks([1,2,3,4,5,6], ['Animal', 'Environmental Sample', 'Human', 'Plants and Products', 'Other', 'Not Provided'])
# Adding  title and labels
plt.title("Number of prophages in bacteria depending on isolation source")
plt.xlabel("Isolation Source")
plt.ylabel("Number of Prophages in Bacteria Genome")

# show plot
plt.show()

