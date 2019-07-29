##Code to create a .csv file to be fed into MYSQL
#the output will be a table of cell, gene and gene count for every cell and every gene
#Input is the merged_out_mmquant.csv file from the PigeonPigs pipeline

output = open('count_for_MySQL.csv',"w")
with open('Merged_out_mmquant.csv',"r") as mmquant:
	gene_list = []
	file = mmquant.readlines()
	cell_list = file[0].split()
	cell_list = cell_list[1:]
	
	#remember which line we are on
	i = 0
	for line in file[1:]:
		gene_list.append(line.split()[0])
		each_line = line.split()[1:]
		
		#remember the position in the list
		j = 0
		for item in each_line:
			if float(item) != 0:
				output.write(cell_list[j] + '\t' + gene_list[i] + '\t' + item + '\n')
			j = j + 1
		i = i + 1
output.close()
