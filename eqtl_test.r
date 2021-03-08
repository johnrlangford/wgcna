library(MatrixEQTL)
exprs <- read.table("BritBams_normalized_counts.txt",header=T,sep="\t",
                    row.names=1)
snps <- read.table("GV_Brit_chr10_GT_Matrix.txt",header=T, sep="\t",
                   row.names=1, nrows=100)


exprs <- exprs[1:100,names(snps)]

pvOutputThreshold = 1e-3
errorCovariance = numeric()
useModel = modelLINEAR

snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns=1
snps$fileSliceSize=2000
snps$LoadFile("GV_Brit_chr22_GT_Matrix.txt")

gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile("BritBams_normalized_counts.txt")

me = Matrix_eQTL_engine(
  snps=snps,
  gene=gene,
  #cvrt=cvrt,
  output_file_name = NULL,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose=TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
);
