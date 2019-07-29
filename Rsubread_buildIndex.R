
#Here is the R install step we used
#BiocManager::install("Rsubread")

#import Rsubread
library(Rsubread)

#Builds the index of human-genome-reference, running Rsubread in build mode
Rsubread::buildindex(basename=file.path("/home/graham/Downloads/our_pipeline/","38_index"), 
                     reference = file.path("/home/graham/Downloads/our_pipeline/GRCh38.primary_assembly.genome.fa"),
                     gappedIndex = TRUE)

#Output files will be prefixed 38_index, allowing targetting by Rsubread Align