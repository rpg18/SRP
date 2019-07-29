###script to find the matching cells in Our_count_matrix (partial output of our attempt to replicate their pipeline) in Their_count_matrix (a downloaded copy of their matrix) and write them to their_67 to allow comparison

#list to hold matching GSM references
GSM_list = []

#list to hold index numbers of matching GSM refs (i.e. their column number)
index_list = []

#Pop open our output to grab the GSM refs from the first line.
with open("Our_count_matrix.csv","r") as Our:
	i = 0
	for line in Our:
		if i>0:
			break
		our_GSM = line.split()
		i = +1
#Open their matrix to check our GSMs against their GSMs
with open("Their_count_matrix.csv", "r") as Their:
	i = 0
	for line in Their:
		their_GSM = line.split()
		j = 0
		for GSM in their_GSM:
			if GSM in our_GSM:
				GSM_list.append(GSM)
				index_list.append(j)
			j = j+1
		if i>0:
			break
		i = +1
#Write matching outputs to file
with open("Their_count_matrix.csv", "r") as Their:
	with open("their_67.csv","w") as Output:
		for line in Their:
			each_line = line.split('\t')
			Output.write(each_line[0]+"\t")
			for index in index_list:
				Output.write(each_line[index]+"\t")
			Output.write("\n")
