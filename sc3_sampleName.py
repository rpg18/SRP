## Generates a new matrix with appropiate sample names for running SC3 with the article's matrix
import os

entries = os.listdir("Data/")

info = [] # initiate info list
genes = [] # initiate gene list
with open('GSM1657871_1772078217.C03.csv') as tsvfile:
    reader = tsvfile.readlines() # read lines from the reads counting file to extract gene names
    for row in reader:
        tri = row.rstrip('\n') # remove new line symbols from data
        th = tri.split('\t') # split data at the tab
        genes.append(th[0]) # append just gene names


for entry in entries: # a for loop to go through list of cells
    with open(entry) as en: # open each data file containing gene counts of each cell
        tr = en.readlines() # read lines of open file
        i=0 # initiate counter
        for rt in tr: # a for loop to go through each line of the of the file
            po = rt.rstrip('\n') # remove new line character
            ug = po.split('\t') # remove tab character
            genes[i] += '\t' + ug[1] # add gene count to correct string
            i+=1 # add one to the counter to move along gene list


# Appends just the sample_name, without chip_information and suffix ".csv" file format 
for entry in entries:
	sample_name = entry.split('_') # split sample name
	info.append(sample_name[0]) # append the prefix: Sample_Name without chip and ".csv" information        


ent = ['Gene_name'] + info # column name in output1.tsv file

with open('output1.tsv', 'w') as ref: # create output1.tsv file
    ref.write('\t'.join(ent)+ '\n') # write columns
    for name in genes: # for loop to write gene count into the file
        ref.writelines(name + '\n')



