''' a code to obtain the information of cell types and append to individual lists and output as one file'''

with open('pipe2_SRR.tsv.csv') as inf: # open information file
	info = inf.readlines() # read information file
	neur=[] # initiate cell type lists
	micr=[]
	end=[]
	fequi=[]
	ferep=[]
    # go through lines in information file
	for row in info: 
		tri = row.rstrip('\n')
		th = tri.split(',')
        # append cell identifiers to cell list depending on cell type
		if 'neurons' in th[2]:
			neur.append(th[0])
		elif 'micro' in th[2]:
			micr.append(th[0])
		elif 'endo' in th[2]:
			end.append(th[0])
		elif 'qui' in th[2]:
			fequi.append(th[0])
		elif 'rep' in th[2]:
			ferep.append(th[0])


with open('raquel.tsv','w') as tr: # open file to write cell lists
    # go through each list of cells and write to list
	tr.write('neurons'+'\t')
	for theo in neur:
		tr.write(str(theo) + '\t')
	tr.write('\n')
	tr.write('microglia'+'\t')
	for qw in  micr:
		tr.write(str(qw) + '\t')
	tr.write('\n')
	tr.write('endothelial'+'\t')
	for er in  end:
		tr.write(str(er) + '\t')
	tr.write('\n')
	tr.write('feotalqui'+'\t')
	for yt in  fequi:
		tr.write(str(yt) + '\t')
	tr.write('\n')
	tr.write('foetalrep'+'\t')
	for ut in ferep:
		tr.write(str(ut) + '\t')
	tr.write('\n')
