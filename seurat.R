library(dplyr) # version 0.8.3
library(Seurat) # version 3.0.2
## Input ----
# Raw reads counting matrix:
trial<- read.csv(("/home/rpg18/Desktop/SRP_article/PIPELINE/Pipeline_SecondSteps/Merged_out_mmquant.csv"), sep='\t')
# Extra article's information:
pipe<- read.csv(("/home/rpg18/Desktop/SRP_article/PIPELINE/Pipeline_SecondSteps/pipe2_SRR.tsv"), sep='\t')

# Matrix preparation for sce object creation
trial_genes <- trial$Gene
rownames(trial) <- trial_genes
trial <- subset(trial, select = -Gene)
dim(trial)

## Initialize the Seurat object with the raw (non-normalized data) ----
pbmc <- CreateSeuratObject(counts = trial, min.cells = 3, min.features = 200)
pbmc
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## Visualize QC metrics as a violin plot ----
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# We decided not to filter our matrix
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Normalising data ----
pbmc <- NormalizeData(pbmc)


## Identification of highly variable features (feature selection) ----
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
CombinePlots(plots = list(plot1, plot2))

## Sacaling the data ----
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## Perform linear dimensional reduction ----
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

## Determine the dimensionality of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20) # drop-off in significance between 9-12

# Alternative heuristic method:
ElbowPlot(pbmc) # the 'elbow' is around 9 and 11 -> the majority of true signal is captured in the first 10 PC's

## Clustering Seurat v3 ----
pbmc <- FindNeighbors(pbmc, dims = 1:10) # 
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 10)