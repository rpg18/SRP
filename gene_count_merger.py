
''' a program to take 466 cells files each containing a column of gene names and a column of gene counts and turning it into a gene count matrix for 466 cells and outputting it into a file'''

import os  #import os module to return the contents of a directory

'''retrieve file names for all cell files'''
entries = os.listdir("Data/")

genes = [] # initiate gene list

with open('Data/GSM1657871_1772078217.C03.csv') as tsvfile: # open file to read genes
    reader = tsvfile.readlines() # read file
    for row in reader: # go through lines in file
        tri = row.rstrip('\n')
        th = tri.split('\t')
        genes.append(th[0]) #append gene names to gene list

for entry in entries: # go through cell files
    with open(entry) as en: # open each cell file
        tr = en.readlines() # read cell file
        i=0 # counter
        for rt in tr: #go through lines of cell file
            po = rt.rstrip('\n')
            ug = po.split('\t')
            genes[i] += '\t' + ug[1] # append counts to genes list
            i+=1

ent = ['Gene_name'] + entries # append column header to the start of cell file list

with open('outputsm.tsv', 'w') as ref: # open file to write to
    ref.write('\t'.join(ent)+ '\n')
    for thin in genes: # for strings in genes list
        ref.writelines(thin + '\n') # write each string to file
