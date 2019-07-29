#
# A code to run the analysis on the dataset produced by Darmanis et al 2015. 
# Data merged and information collated from gene_cell_merger.py
# including tsne, minimum spanning tree, pca, MHC boxplots and top 20 genes by cluster
#

# load packages
library(scde)
library(mclust)
library(tsne)
library(parallel)
library(igraph)
library(FactoMineR)
library(tibble)
library(plotly)
library(ggplot2)


# read in tsv files created by gene_merge_count.py
wholematrix <- read.delim('~/SteeredRP/wholematrix.tsv', header=TRUE, row.names = 1)
neurons <- read.delim("~/SteeredRP/neurons.tsv")
foetal <- read.delim("~/SteeredRP/prenatal.tsv")
endothelial<-read.delim("~/SteeredRP/endothelial.tsv")
microglia<-read.delim("~/SteeredRP/microglia.tsv")
fetalqui<- read.delim("~/SteeredRP/fetalqui.tsv")
fetalrep<-read.delim("~/SteeredRP/fetalrep.tsv")

# create matrix of specified cells
neuronmat<-wholematrix[,neurons$V1] 
foetalmat<-wholematrix[,foetal$V1]

# create a matrix of fetal and adult neurons
fenemerge <- merge(foetalmat, neuromat, by=0, all = TRUE)
row.names(fenemerge)<-fenemerge$Row.names
fenemerge$Row.names<-NULL

#a function to sort the matrix, removing cells below threshold, creating a distance matrix of the cells
sorting<-function(fileread){
  cd<-clean.counts(fileread,min.lib.size=1800,min.reads=1,min.detected=1)
  o.ifm <- scde.error.models(counts = cd, n.cores = 2, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
  valid.cells <- o.ifm$corr.a > 0
  o.ifm <- o.ifm[valid.cells, ]
  o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
  
  p.self.fail <- scde.failure.probability(models = o.ifm, counts = cd)
  n.simulations <- 1000; k <- 0.9
  cell.names<-colnames(cd)
  names(cell.names) <- cell.names
  dl <- mclapply(1:n.simulations,function(i) {
    scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
      x <- cd[,nam];
      x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
      x;
    }))
    rownames(scd1) <- rownames(cd); 
    cor(log10(scd1+1),use="pairwise.complete.obs");
  }, mc.cores = 3)
  there<-1-Reduce("+",dl)/length(dl)
  return(there)}

# a function to run tsne and MClust on a distance matrix and plot BIC graph as well as the accompanying boxplot, plot a 3D tsne coloured by cluster and cell type
tsned<-function(distmat){
  direct.dist<- as.dist(distmat)
  sne <-tsne(direct.dist, k=3)
  yoghurt<-Mclust(sne, G=1:40, modelNames = c("EII","VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV"),prior=priorControl())
  plot(yoghurt$BIC, colors = c("gray", "black", "orange", "darkred", "red", "magenta", "darkgreen", "green","lightblue","darkblue"))
  boxplot(yoghurt$z)
  cluslist<- c()
  for (num in yoghurt[["classification"]]){
    if (num %in% 1){cluslist<-append(cluslist,'black')}
    else if (num %in% 2){cluslist<-append(cluslist,'orange')}
    else if (num %in% 3){cluslist<-append(cluslist,'purple')}
    else if (num %in% 4){cluslist<-append(cluslist,'green')}
    else if (num %in% 5){cluslist<-append(cluslist,'red')}
    else if (num %in% 6){cluslist<-append(cluslist,'blue')}
    else if (num %in% 7){cluslist<-append(cluslist,'yellow')}
    else if (num %in% 8){cluslist<-append(cluslist,'cyan')}
    else if (num %in% 9){cluslist<-append(cluslist,'pink')}
  }
  grahamlist<-c()
  for(k in row.names(all)){
    if (k %in% opc$V1){grahamlist<-append(grahamlist, 'green')}
    else if(k %in% astrocyte$V1){grahamlist<-append(grahamlist, 'purple')}
    else if(k %in% micro$V1){grahamlist<-append(grahamlist, 'blue')}
    else if(k %in% endo$V1){grahamlist<-append(grahamlist, 'black')}
    else if(k %in% neuls$V1){grahamlist<-append(grahamlist, 'orange')}
    else if(k %in% fetalqui$V1){grahamlist<-append(grahamlist, 'red')}
    else if(k %in% replfe$V1){grahamlist<-append(grahamlist, 'magenta')}
    else if(k %in% hybrid$V1){grahamlist<-append(grahamlist, 'cyan')}
    else if(k %in% oligo$V1){grahamlist<-append(grahamlist, 'grey')}
    else{grahamlist<-append(grahamlist, k)}
  }  
  
  plot_ly(as.data.frame(sne), x = sne[,1], y= sne[,2], z= sne[,3], type = 'scatter3d', marker= list(color = cluslist, size=1))
  plot_ly(as.data.frame(sne), x = sne[,1], y= sne[,2], z= sne[,3], type = 'scatter3d', marker= list(color = grahamlist, size=1))
}

# a function to create the minimum spanning tree object with colour by community
tree<-function(distmat){
  grapho<-graph.adjacency(as.matrix(distmat), weighted=TRUE, mode = "lower")
  minm<-mst(grapho)
  com<-cluster_walktrap(minm)
  V(minm)$color<-c("blue", "green", "magenta", "red", "gray", "black", "orange", "darkred", "magenta", "darkgreen","lightblue","darkblue", "purple")[membership(com)] #add 13 cols to list
  return(minm)
}

# a function to perform PCA and assign colour by previously assigned cell type and plot
pca<-function(distmat, file1, file2, file3){
  gre<-PCA(distmat, graph = FALSE) # make pca not plot random plots
  collist = c()
  for(k in row.names(distmat)){
    if (k %in% file1[,1]){
      print(k)
      collist<-append(collist, 'green')
    }
    else if (k %in% file2[,1]){
      collist <- append(collist, 'purple')
    }
    else if (k %in% file3[,1]){
      collist <- append(collist, 'red')
    }
  }
  plot_ly(as.data.frame(gre$ind$coord), x = ~Dim.1, y= ~Dim.2, z= ~Dim.3,type = 'scatter3d', marker= list(color = collist, size=1))
  return(gre)
}

#a function to log counts per million normalise data
norm <-function(data){
  for (i in 1:length(names(data))){
    sumf<-sum(data[,i])
    for (cell in 1:length(data[,i])){
      cell1 = (data[cell,i]/sumf)*10^6
      if (cell1 > 1) {
        cell2 = log10(cell1)
      }
      else{
        cell2 = log10(1)
      }
      data[cell,i] <- cell2
    }}
  return(data)}

# a function to create smaller matrix with only specified genes and cells
matrix<-function(input, matrix){
  genesss <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAPBP", "CALR", "ERAP1", "HSPA5", "PDIA3", "TAP2", "SEC61A1", "SEC61A2", "SEC61B", "SEC61G")
  cols<-as.vector(input$V1)
  matrix1<-matrix[genesss,input]
  return(matrix1)
}

#create distance matrices for different collection of cell
allsort<-sorting(wholematrix)
neuronsort<-sorting(neuronmat)
foetalsort<-sorting(foetalmat)
mergesort<-sorting(fenemerge)

# run tsne on entire distance matrix
altsne<-tsned(allsort)

# create fetal and adult neuron minimum spanning trees
fetreegra<-tree(foetalsort)
nertreegra<-tree(neuronsort)

# plot adult and fetal minimum spanning tree
plot(fetreegra, layout=layout.fruchterman.reingold, vertex.size=4, vertex.label=NA, asp=FALSE, edge.arrow.mode=0)
plot(nertreegra, layout=layout.fruchterman.reingold, vertex.size=4, vertex.label=NA, asp=FALSE, edge.arrow.mode=0)

#run pca fucntion on combined adult and fetal neurons
fenepca<-pca(mergesort, quiefe, replfe, neuls)

# Code to plot MHC boxplots 

#normalise whole matrix
normwholematrix<-norm(wholematrix)

#make smaller matrices from selected cells
gennormendomat<-matrix(endothelial, normwholematrix)
gennormfequimat<-matrix(fetalqui, normwholematrix)
gennormferepmat<-matrix(fetalrep, normwholematrix)
gennormmicromat<-matrix(microglia, normwholematrix)
gennormneuromat<-matrix(neurons, normwholematrix)

# identify gene names
name <- rownames(genenormneuromat)

#convert matrices to dataframe and add colnames (genes)
data1 <- as.data.frame(t(genenormneuromat))
colnames(data1) <- name
data2 <- as.data.frame(t(gennormendomat))
colnames(data2) <- name
data3 <- as.data.frame(t(gennormmicromat))
colnames(data3) <- name
data4 <- as.data.frame(t(gennormfequimat))
colnames(data4) <- name
data5 <- as.data.frame(t(gennormferepmat))
colnames(data5) <- name

# make one large dataset
All <- rbind(data1, data2, data3, data4, data5)
Sample <- c(rep("Adult neurons", nrow(data1)),rep("Endothelial cells", nrow(data2)),rep("Microglia", nrow(data3)),
            rep("Fetal neurons (quiescent)", nrow(data4)),rep("Fetal neurons (replicating)", nrow(data5)))

#identify groups and attach data
All <- add_column(All, Sample, .before= All$`HLA-A `)
attach(All)

#identify factors to group data
All$Sample<-factor(All$Sample, levels = c("Adult neurons","Endothelial cells","Microglia","Fetal neurons (quiescent)","Fetal neurons (replicating)"))
All$Sample <- as.factor(All$Sample)

# calculate mean of SEC61 cells
All$SEC61<-rowMeans(All[,12:15])

#plot boxplot of each gene
HLA_Box <- ggplot(data=All, aes(x=Sample, y=`HLA-A `)) + geom_boxplot()
HLA_Box + geom_jitter(shape = 19, position = position_jitter(0.3))+ labs(title="HLA-A",x="Type of cell",y="LogCPM")

#Code to find the top 20 genes

#initiate cluster lists
cluster1 <-c()
cluster2 <-c()
cluster3 <-c()
cluster4 <-c()
cluster5 <-c()
cluster6 <-c()
cluster7 <-c()
cluster8 <-c()
cluster9 <-c()

# append cells to cluster lists
count = 1
for(k in row.names(all)){
  meh <- yoghurt$classification[count]
  count <- count + 1
  if (meh %in% 1){cluster1<-append(cluster1,k)}
  else if (meh %in% 2){cluster2<-append(cluster2,k)}
  else if (meh %in% 3){cluster3<-append(cluster3,k)}
  else if (meh %in% 4){cluster4<-append(cluster4,k)}
  else if (meh %in% 5){cluster5<-append(cluster5,k)}
  else if (meh %in% 6){cluster6<-append(cluster6,k)}
  else if (meh %in% 7){cluster7<-append(cluster7,k)}
  else if (meh %in% 8){cluster8<-append(cluster8,k)}
  else if (meh %in% 9){cluster9<-append(cluster9,k)}
}
#create a gene count matrix for each cluster
hopeful1<-wholematrix[,cluster1]
hopeful2<-wholematrix[,cluster2]
hopeful3<-wholematrix[,cluster3]
hopeful4<-wholematrix[,cluster4]
hopeful5<-wholematrix[,cluster5]
hopeful6<-wholematrix[,cluster6]
hopeful7<-wholematrix[,cluster7]
hopeful8<-wholematrix[,cluster8]
hopeful9<-wholematrix[,cluster9]

#a function to identify the top 20 genes in each cluster
topgene <- function(hopeful){
  total<- rowSums(hopeful);
  hopeful<-cbind(hopeful,total);
  ordered1 <- hopeful[order(-total),];
  topgenes <-row.names(ordered1[1:20,]);
  return(topgenes)
}

#finding the top 20 genes for each cluster
genes1<-topgene(hopeful1)
genes2<-topgene(hopeful2)
genes3<-topgene(hopeful3)
genes4<-topgene(hopeful4)
genes5<-topgene(hopeful5)
genes6<-topgene(hopeful6)
genes7<-topgene(hopeful7)
genes8<-topgene(hopeful8)
genes9<-topgene(hopeful9)

#top 20 genes overall
out<-topgene(wholematrix)

