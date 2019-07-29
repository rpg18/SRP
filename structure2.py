''' a program to take 466 cells files each containing a column of gene names and a column of gene counts and turning it into a gene count matrix for 466 cells and outputting it into a file, and cells list identifying different cell types and outputting each to a file, adult neuronal, endothelial, microglia, oligodendrocytes, astrocytes, OPC, hybrid, prenatal, fetal neurons quiescent and fetal neurons replicating '''

import os # import os module to return the contents of a directory

''' a function to obtain the list of genes that will for the row names of the matrices, output gene list'''
def geneList():
    genes = [] # initiate gene list
    with open('Data/GSM1657871_1772078217.C03.csv') as tsvfile: # open one data file in directory
        reader = tsvfile.readlines() # read each line of the data file
        for row in reader: # a for loop to go through each row of the data file
            tri = row.rstrip('\n') # remove new line symbols from data
            th = tri.split('\t') # split data at the tab
            genes.append(th[0]) #append gene names to gene list
    return(genes) # return gene list

''' a function to get the gene count of each gene from each cell and append to list, input 'entries' a list of cells, 'genesf' a list of genes, output a list of strings of the counts to each gene'''
def entryList(entries, genesf):
    for entry in entries: # a for loop to go through list of cells
        with open('Data/' + entry) as en: # open each data file containing gene counts of each cell
            tr = en.readlines() # read lines of open file
            i=0 # initiate counter
            for rt in tr: # a for loop to go through each line of the of the file
                po = rt.rstrip('\n') # remove new line character
                ug = po.split('\t') # remove tab character
                genesf[i] += '\t' + ug[1] # add gene count to correct string
                i+=1 # add one to the counter to move along gene list
    return(genesf) # return list of strings

''' a function to trim the last three rows from the the matrix (no_feature, ambiguous, alignment_not_unique)'''
def trim(fre):
    for i in range(3): # a loop of three
        fre.pop(-1) # remove row
    return fre # return new matrix

''' a function to retrieve information to label cells'''
def runTable():
    thing = [] # initiate list to fill with cell information
    with open('SraRunTable.tsv') as filed: # open SraRunTable file containing information
        readr = filed.readlines() # read each line of the file
        for ro in readr: # for to go through every line of the file 
            number = ro.split('\t') # split each row with a tab
            thing.append(number[9:13]) # append information to list
    return thing #return information list

''' a function to create lists to sepate foetal and neuronal cells to create separate matrices'''
def seParate(thing, entry):
    # initiate names lsit for each cell type
    prenat = [] 
    neur=[] 
    fequi=[]  
    ferep=[] 
    micr=[]
    end=[]
    olig=[]
    hybrid=[]
    astro=[]
    opc=[]
    for entry in entries: # for loop through list of cells
        for item in thing: # for loop through information list
            if item[0] in entry: # if statement to compare cell identifiers
                #append cell identify to cell type list accordingly
                if 'prenatal' in item[1]:
                    prenat.append(entry)
                if 'neurons' in item[3]:
                    neur.append(entry)
                elif 'micro' in item[3]:
                    micr.append(entry)
                elif 'endo' in item[3]:
                    end.append(entry)
                elif 'qui' in item[3]:
                    fequi.append(entry)
                elif 'rep' in item[3]:
                    ferep.append(entry)
                elif 'olig' in item[3]:
                    olig.append(entry)
                elif 'hybrid' in item[3]:
                    hybrid.append(entry)
                elif 'astro' in item[3]:
                    astro.append(entry)
                elif 'OPC' in item[3]:
                    opc.append(entry)
    return prenat, neur, fequi, ferep, micr, end, olig, hybrid, astro, opc

''' a function that takes a filename which data will be outputted into, a list of genes and a list of cells. Outputting a matrix containing cells and their cells counts for all of the genes into given file.'''

def outPut(filename, genesk, names):
    with open(filename, 'w') as ref:
        ref.write('\t'.join(names)+ '\n')
        for thin in genesk:
            ref.writelines(thin + '\n')


'''A function for to return a file with a list of cells in, one to a line. Takes a string for a filename and list of cell names. '''

def fileinp(filename, lists):
    with open(filename, 'w') as fr:
        for naomi in lists:
            fr.writelines(str(naomi) + '\n')

'''retrieve file names for all cell files'''
entries = os.listdir("Data/")
entries.sort() #put entries in alphanumeric order
names = [''] + entries

'''create a list of genes'''
genes = geneList()

'''create a list of strings '''
entry = entryList(entries, genes)

trimmed = trim(entry) #trim unecessary information from list of strings

table = runTable() #extra information from supplementary table

prenata, neuro, fequil, ferepl, micro, endo, oligo, hybrid, astro, opc = seParate(table, entries) # separate cell names

outPut('wholematrix.tsv', trimmed, names) # output cell count matrix to file

#output lists of cell types to files
fileinp('prenatal.tsv', prenata)
fileinp('neurons.tsv', neuro)
fileinp('fetalqui.tsv', fequil)
fileinp('fetalrep.tsv', ferepl)
fileinp('microglia.tsv', micro)
fileinp('endothelial.tsv', endo)
fileinp('oligo.tsv', oligo)
fileinp('hybrid.tsv', hybrid)
fileinp('astrocyte.tsv', astro)
fileinp('opc.tsv', opc)

