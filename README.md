# SRP


# Re-Analysis pipeline: ra
1- sra_download.pl \
2- first_untrimmed_fastq.sh \
3- STAR_build.sh \
4- second_trimmed_fastq.sh \
5- structure2.py\
6- original_data_analysis.R

# PigeonPigs pipeline: pp 
1- Rsubread_bioldIndex \
2- Rsubread.R \
2.1- Rsubread_mmquant.sh  
3- sc3_sce_Info.py \
3.1- sc3_inputDarmanis.py \
4- sc3_clustering.R \
5- seurat.R

# Website: wb
1- FOLDER: Wobsite \
2- FOLDERS: BIC, Diffusion_Map, MST, PCA_neurons_by_cell_type, gene_boxplot, our_tSNE_PCA_sce, plotDcater, tSNE \
3- P2XRbase.php \
4- project_sql.sql 

