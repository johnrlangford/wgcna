setwd("/Users/john/projects/wgcna")
library(HelpersMG)
library(WGCNA)
library(tidyverse)

#important setting according to documentation
options(StringsAsFactors = FALSE)

expr = read.table("BritBams_normalized_counts.txt", row.names=1)
expr = t(expr)
expr = log2(expr+1)
geneMeans = apply(as.matrix(expr),FUN=mean, MARGIN=2)

#log of expression seems to give better results
#not necessary though

#get n most highly expressed genes
n=1000
highexprs = tail(sort(geneMeans),n=n)

expr = expr[ , names(highexprs)]

#cluster and look for outliers

sampleTree = hclust(dist(expr), method="average")

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#dendrogram shows 1 outlier HG00102

clust = cutree(sampleTree, h=22)
table(clust)

keepSamples = clust==1

expr = expr[keepSamples, ]

nGenes = ncol(expr)
nSamples = nrow(expr)

save(expr, file="BritBams_cleaned_counts.RData")

#next line turns on multithreading. doesn't work in rstudio
#enableWGCNAThreads()

#raise coexpression similarity to power to calculate adjacency
#following picks power
powers = c(c(1:10), seq(from=12, to=20, by=2))
sft = pickSoftThreshold(expr, powerVector=powers, verbose=5)
par(mfrow=c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#pick first power that gives approximate scale independence
#net is reference to module object
net = blockwiseModules(expr, power = 6,
                       TOMType = "unsigned", minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "britBamTOM",
                       verbose = 3)

#examine net object. colors correspond to modules. 0 is for unassigned genes
table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs,moduleLabels,moduleColors,geneTree,
     file="britBamNetworkConstruction_auto.RData")
