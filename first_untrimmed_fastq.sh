##PIPELINE for processing Raw RNA data as processed in paper

#Requires the fafapa_call.py and fafapa_callCAPS.py scripts in the same working directory

#Code for 1st STAR run to generate the indexed version of genome, needs to be run on something with lots of RAM and speed: STAR --runThreadN 12 --runMode genomeGenerate --genomeDir hg19_index --genomeFastaFiles hg19.fa --sjdbGTFfile hg19.GTF --sjdbOverhang 74 --genomeSAsparseD 4

#prefix for the path of the Raw RNA Data directory (containing subdirectories, which each contain 2 fastq.gz files (the paired reads for each cell)
#Enter the FULL path here to prevent STAR throwing an error
path="/home/graham/Downloads/trial/Data/"
#loop through subdirectories in Directory 
for i in $path*
do
	#suffixes for paired reads
	one="_1.fastq.gz"
	two="_2.fastq.gz"
	#j becomes the current file path, minus the dir path - i.e. the name of the subdir
	j=${i#$path}
	#one and two are concatenated from previous parts to create the full path to the fastq files
	one=$i"/"$j$one
	two=$i"/"$j$two

	#runs prinseq++ as specified in paper on the 2 fastq files
	prinseq++ -min_len 30 -trim_left 10 -trim_qual_right 25 -lc_entropy 0.65 -threads 1 -fastq $one -fastq2 $two -out_name $j -out_bad /dev/null/ -out_bad2 /dev/null/ -out_single /dev/null/ -out_single2 /dev/null/
	
	#set the suffixes to the new output pair
	one="_good_out_R1.fastq"
	two="_good_out_R2.fastq"

	#runs FASTQC as specified in paper on prinseq output
	fastqc $j$one $j$two --extract -t 1

	#calls a python script that parses FASTQC output to find overrepresented sequences for cutadapt (first pair)
	parse_fastqc1=$(python3 fafapa_call.py $j"_good_out_R1")

	#calls a python script that parses FASTQC output to find overrepresented sequences for cutadapt (2nd pair)
	parse_fastqc2=$(python3 fafapa_callCAPS.py $j"_good_out_R2")

	#runs cutadapt as specified in paper, removing overrepresented sequences from FASTQC from prinseq output
	cutadapt $parse_fastqc1 $parse_fastqc2 -e 0.15 -m 30 -o cutadR1$j.fastq -p cutadR2$j.fastq  -j 1 $j$one $j$two

	#runs prinseq++ as specified in paper on cutadapt output
	prinseq++ -min_len 30 -fastq cutadR1$j.fastq -fastq2 cutadR2$j.fastq -out_name trimmed$j -out_bad /dev/null/ -out_bad2 /dev/null/ -out_single /dev/null/ -out_single2 /dev/null/

	#runs Trim Galore as specified in paper on 2nd prinseq output
	trim_galore trimmed$j$one trimmed$j$two --stringency 1 --nextera --paired

	#Feed the trimmed files into STAR (with settings as in paper. runThreadN 4 to use all 4 cores on laptop.
	#STAR likes full file paths, so manual entry of genomeDir (where the aligned genome is) is required
	STAR --readFilesIn "trimmed"$j"_good_out_R1_val_1.fq" "trimmed"$j"_good_out_R2_val_2.fq" --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMstrandField intronMotif --runThreadN 3 --genomeDir "/home/graham/Downloads/trial/hg19_index" --outFileNamePrefix $j

	#Create a gene count matrix for the current cell with htseq (setting as in paper) and write to output
	htseq-count -m intersection-nonempty -s no -f sam $j"Aligned.out.sam" UCSC_geneid.gtf > $j"_out.csv"

	#These can now be merged with the gene_count_merger.py script
done
