###Parses fastqc output to find overrepresented sequences to feed into cutadapt

# imports
#may need to install Fadapa (pip install fadapa)
from fadapa import Fadapa
import sys

#take argument from bash script (which will be $j - the UID of the cell)
name = sys.argv[1]

#load file into fadapa parser
f = Fadapa('/home/graham/Downloads/trial/'+name+'_fastqc/fastqc_data.txt')

#get raw data for Overrepresented sequences
pass_seq = f.raw_data('Overrepresented sequences')[0]

#Initialise list of seqs
list_of_seqs = []

#If there are no overrepresented sequences, the clean parser breaks!
#Therefore, we cannot reference .clean unless .raw contains something
if pass_seq != ">>Overrepresented sequences	pass":
	#Loop through the .clean parsed data
	for data in f.clean_data('Overrepresented sequences'):
		#Add the first index of the .clean data to list
		#First entry will by #Sequence, subsequent will be the actual seqs
		list_of_seqs.append(data[0])

#Create empty output string
output = ""

#Loop through the list of sequeces from index 1 onwards (as the index 0 will be #Sequence)
for things in list_of_seqs[1:]:
	#Concatenates seqs in format required for cutadapt argument (first pair)
	output = output + " -a " + things

#Prints the cutadapt argument - this sends it the stdout, allowing it to feed into the bash script
print(output))
