## This script creates two tsv files that store sample name (SRR or GSM prefixes)
# pipe2_GSM is the file that we need to generate the SingleCellExperiment (SCE) object: For Darmaris et al, 2015 Matrix
# pipe2_SRR.tsv is the file that we need to generate the SCE object: For MMquant Matrix 

liston =[] # open liston list
with open("SraRunTable.txt") as file: # data from Darmanis et at, 2015
    lines = file.readlines()	# read lines from the text file
    for row in lines:	# a for loop to go through rows of the file
        data = row.split('\t') # split data at the tab
        liston.append(data) # append data into liston list
		
        with open('pipe2_SRR.tsv', 'w') as ref: # create and write new file
            info = [[aa[7],aa[10],aa[12]] for aa in liston]
			# list of [SRR_sample_name, age of individuals and cell_type]
            for item in info:	# a for loop to go through items of info list              
                dataSRR = '\t'.join(item) # join strings at the tab
                ref.writelines(dataSRR+'\n') # write strings-rows into new file


        with open('pipe2_GSM.tsv', 'w') as ref: # create and write new file
            info = [[aa[9],aa[10],aa[12]] for aa in liston]
			# list of [GSM_sample_name, age of individuals and cell_type]
            for item in info:	# a for loop to go through items of info list            
                dataGSM = '\t'.join(item) # join strings at the tab
                ref.writelines(dataGSM+'\n') # write strings-rows into new file



