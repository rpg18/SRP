#Code for 1st STAR run to generate the indexed version of genome 
#needs to be run on something with lots of RAM and CPU
#Requires an output directory (hg19_index), and inputs of a genome assembly (hg19.fa) and genome annotation file (hg19.GTF)

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir hg19_index --genomeFastaFiles hg19.fa --sjdbGTFfile hg19.GTF --sjdbOverhang 74 --genomeSAsparseD 4
