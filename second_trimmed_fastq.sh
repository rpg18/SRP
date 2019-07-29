##PIPELINE for processing Raw RNA data as processed in paper

#Requires a Data directory containing subdirs named for each cell, an hg19_index directory containing the output of the genome align STAR run, and a UCSC_geneid.gtf file containing the proper gene names for each gene

#Code for 1st STAR run to generate the indexed version of genome, needs to be run on something with lots of RAM and speed: STAR --runThreadN 12 --runMode genomeGenerate --genomeDir hg19_index --genomeFastaFiles hg19.fa --sjdbGTFfile hg19.GTF --sjdbOverhang 74 --genomeSAsparseD 4

#Outputs logs and a sam file for the last cell processed (irrelevant), and a gene expression matrix for each cell processed

#prefix for the path of the Raw RNA Data directory (containing subdirectories, which each contain 2 fastq.gz files (the paired reads for each cell)
#Enter the FULL path here to prevent STAR throwing an error
path="/home/graham/Downloads/trial/Data/"
#loop through files in Dir
for i in $path*
do
	#suffixes for paired reads in their gzip files
	one="_1.fastq.gz"
	two="_2.fastq.gz"

	#j becomes the current file path, minus the dir path - i.e. the name of the subdir
	j=${i#$path}

	#one and two are concatenated from previous parts to create the full path to the fastq files
	one=$i"/"$j$one
	two=$i"/"$j$two

	#Feed the files into STAR (with settings as in paper. runThreadN 4 to use all 4 cores on laptop.
	#STAR likes full file paths, so manual entry of genomeDir (where the aligned genome is) is required
	STAR --readFilesIn $one $two --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMstrandField intronMotif --runThreadN 4 --genomeDir "/home/graham/Downloads/trial/hg19_index" --readFilesCommand zcat

	#Feed the files into htseq to produce expression per gene for each cell
	htseq-count -m intersection-nonempty -s no -f sam "Aligned.out.sam" UCSC_geneid.gtf > $j"_out.csv"

	#These can now be merged with the gene_count_merger.py script
done
