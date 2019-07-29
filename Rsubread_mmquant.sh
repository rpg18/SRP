###Code to produce gene count matrix for our pipeline
#May require an external hard drive to hold the BAM output, as the fastq files take a substantial amount of storage space
#Run in a directory with the Rscript Rsubread 

#Requires a previous run of Rsubread_buildIndex.R script

#prefix for the path of the directory that contains the correct gtf file, the output of an Rsubread_buildIndex.R run (here, the 38_index files), and a Raw_RNA_Data directory containing subdirectories, which each contain 2 fastq.gz files (the paired reads for each cell)
#
path="/home/graham/Downloads/our_pipeline/"

#loop through subdirectories in Raw_RNA_Data Directory
for i in $path"Raw_RNA_Data/"*
do
	#suffixes for paired reads in their gzip files
	fastq1="_1.fastq.gz"
	fastq2="_2.fastq.gz"

	#j becomes the current file path, minus the dir path - i.e. the name of the subdir
	j=${i#$path"Raw_RNA_Data/"}

	#one and two are concatenated from previous parts to create the full path to the fastq files
	one=$i"/"$j$fastq1
	two=$i"/"$j$fastq2

	#calls Rsubread R script on the pairs of fastq files (this runs the align step)
	Rscript Rsubread.R $one $two $path"38_index" $j".bam"

	#moves the BAM file output of the Rscript to an external hard drive
	mv $j".bam" "/media/graham/External2TB/Graham/"

	#concatenate the path to the BAM file onto a string.
	#Finally, this string will hold all file paths for all outputted BAM files
	total=$total" ""/media/graham/External2TB/Graham/"$j".bam"

done

#Runs mmquant on all of the BAM files to generate a count matrix for all genes and cells
# -a holds the path to the gtf file used, -r the string of all BAM file paths from above, -t controls number of threads to use,-g ensures output uses gene names, -o name of gene count matrix ouput file, and -O name of the file contains stats on aligned reads, etc
mmquant -a $path"gencode.v31.annotation.gtf" -r $total -t 4 -g -o "All_Out_mmquant.csv" -O "Stats_report"
echo "Job done!"
#The mmquant count matrix will contain the useful merged genes. It will also contain reads that mapped uniquely to parts of the same gene, and are stored (thanks to the -g option) under the same name on different rows. These counts can now be combined by the merge_dupl_gene.py script.
