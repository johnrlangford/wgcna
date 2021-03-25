library(MatrixEQTL)
SNP_file_name = "GV_Brit_chr19_GT_Matrix.txt"

loadSNPS <- function(SNP_file_name,MEs) {
  snps = read.table(SNP_file_name, row.names = 1, header = T)
  common_samples = intersect(names(MEs), names(snps))
  snps = snps[ , common_samples]
  SNPS = SlicedData$new()
  SNPS$CreateFromMatrix(as.matrix(snps))
}
useModel = modelLINEAR

#no covariates so set covariates filename to character()
covariates_file_name = character()
output_file_name = tempfile()
errorCovariance = numeric()
#instead of modeling gene expression, we model eigengene expression for
#module eigengenes
load("britBamNetworkConstruction_auto.RData")
expression_file_name = MEs
MEs = t(MEs)
MEs = data.frame(MEs)

pvOutputThreshold = 1e-5

SNPS <- loadSNPS(SNP_file_name, MEs)

#reduce snps by sampling to speed up tests

#sample_snps = sample(row.names(snps), size=)

#snps = snps[sample_snps, ]


MEs = MEs[ , common_samples]



eigengenes = SlicedData$new()
eigengenes$CreateFromMatrix(as.matrix(MEs))

#iterate through all snp files
#in directory
#snp file names begin with GV_

directory_files = list.files()
snp_files = directory_files[str_detect(directory_files, regex("^GV"))]

for (val in snp_files) {
  SNPS = loadSNPS(val,MEs)
  chromosome = str_extract(val, "chr\\d+")
  output_file_name = paste(chromosome,"eqtl",sep="_")
  print(output_file_name)
  me = Matrix_eQTL_engine(
    snps=SNPS,
    gene=eigengenes,
    output_file_name=output_file_name,
    pvOutputThreshold=pvOutputThreshold,
    useModel=useModel,
    errorCovariance=errorCovariance,
    verbose=TRUE,
    pvalue.hist=TRUE,
    min.pv.by.genesnp=FALSE,
    noFDRsaveMemory=FALSE
  )
}
me = Matrix_eQTL_engine(
  snps=SNPS,
  gene=eigengenes,
  output_file_name=output_file_name,
  pvOutputThreshold=pvOutputThreshold,
  useModel=useModel,
  errorCovariance=errorCovariance,
  verbose=TRUE,
  pvalue.hist=TRUE,
  min.pv.by.genesnp=FALSE,
  noFDRsaveMemory=FALSE
)

