library(WGCNA)
library(DESeq2)
library(org.Hs.egENSEMBL2EG)
options(stringsAsFactors = FALSE)
expr <- read.table("BritBams_normalized_counts.txt", header = TRUE, sep="\t",
                 row.names=1)
gsg <- goodSamplesGenes(expr)

gsg$allOK


expr <- t(expr)
expr <- data.frame(expr)

#pick random subset of 1000 genes

genes <- sample(names(expr),30000)

expr <- expr[ , genes]
expr <- log2(expr+1)
powers = c(c(1:10), seq(from=12, to=20, by=2))
sft = pickSoftThreshold(expr, powerVector = powers, verbose=5)

par(mfrow = c(1,2))
cex1=.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

cor <- WGCNA::cor
net = blockwiseModules(expr, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
