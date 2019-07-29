# pipe2.tsv is the file that we need to generate the "sce" object, which contains two columns: age and cell_type. For Darmaris et al, 2015 Matrix
# pipe2_SRR.tsv is the file that we need to generate the "sce" object, which contains two columns: age and cell_type. For MMquant Matrix

list =[] # stores all lines from the SraRunTable file

# read SraRunTable where is stored the sample data
with open("SraRunTable.txt") as file:
    lines = file.readlines()
    for row in lines:
        line = row.split('\t')
        list.append(line)

		# pipe2.tsv where we will store the additional information for the creation of the SCE object
        with open('pipe2.tsv', 'w') as ref:
       		info = [[aa[8],aa[9],aa[12]] for aa in list] # stores: Sample_Name (GSM_name), age, cell_type
       		for item in info:
       			feature = '\t'.join(item)
       			ref.writelines(feature+'\n')

		# pipe2_SRR.tsv where we will store the additional information for the creation of the SCE object
        with open('pipe2_SRR.tsv', 'w') as ref:
       		info = [[aa[7],aa[9],aa[12]] for aa in list] # stores: Run (SRR_name), age, cell_type 
       		for item in info:
       			feature = '\t'.join(item)
       			ref.writelines(feature+'\n')
