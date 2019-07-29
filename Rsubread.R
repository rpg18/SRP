#R commands to be called from bash script using Rscript command

#takes args from bash script input
args = commandArgs(trailingOnly=TRUE)

# test if arguments are OK
if (length(args)!=4) {
  stop("Eejit. Needs read1, read2, gtf and output paths.n", call.=FALSE)
}

#the 1st arg is the first fastq file of the pair
fq_R1 = file.path(args[1])
#2nd is the 2nd fastq file of the pair
fq_R2 = file.path(args[2])
#3rd is the output of the Rsubread-Build run
IndexPath = file.path(args[3])
#4th is the output file name
Outpath = file.path(args[4])

#Rsubread run in align mode on the paired files
Rsubread::align(index=IndexPath,
                readfile1=fq_R1,
                readfile2=fq_R2,
                type=0, #0 for rna
                output_file=Outpath,
                nthreads = 4 #Number of cores to use...
)
#This will output a BAM file, which will be fed into mmquant by the bash script
