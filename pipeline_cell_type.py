## This script creates two tsv files that store sample name (SRR or GSM prefixes)

liston =[] # open liston list
with open("SraRunTable.txt") as file: # data from Darmanis et at, 2015
    lines = file.readlines()	# read lines from the text file
    for row in lines:	# a for loop to go through rows of the file
        data = row.split('\t') # split data at the tab
        liston.append(data) # append data into liston list
		
        with open('pipe2_SRR.tsv', 'w') as ref: # create and write new file
            info = [[aa[7],aa[10],aa[12]] for aa in liston]
			# list of [SRR_sample_name, age of individuals and cell_type]
            for item in info:                 
                wow = '\t'.join(item)
                ref.writelines(wow+'\n')


        with open('pipe2.tsv', 'w') as ref: # create and write new file
            info = [[aa[9],aa[10],aa[12]] for aa in liston]
			# list of [GSM_sample_name, age of individuals and cell_type]
            for item in info:                 
                wow = '\t'.join(item)
                ref.writelines(wow+'\n')



