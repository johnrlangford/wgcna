library(MatrixEQTL)

snpDF = read.table("GV_Brit_chr10_GT_Matrix.txt", sep="\t", header = T,
                   row.names=1, nrows=100)
exprDF = read.table("BritBams_normalized_counts.txt", row.names = 1, header=T)
exprDF = exprDF[names(snpDF)]

write.table(exprDF, "expression.txt", sep="\t",row.names=T, col.names=T)
snpFile = "GV_Brit_chr10_GT_Matrix.txt"
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile(snpFile)

#generate small snp file for testing

snps2 = SlicedData$new()

exprFile = "expression.txt"
genes = SlicedData$new()
genes$fileDelimiter = "\t"
genes$fileOmitCharacters = "NA"
genes$fileSkipRows = 1
genes$fileSkipColumns = 1
genes$fileSliceSize = 2000
genes$LoadFile(exprFile)

output_file_name = tempfile()
pvOutputThreshold = 1e-2
useModel = modelLINEAR
errorCovariance = numeric()
covariatesFile = character()
cvrt = SlicedData$new()
cvrt$fileDelimiter = '\t'
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSliceSize = 2000
cvrt$LoadFile(covariatesFile)
me = Matrix_eQTL_engine(
  snps = snps,
  gene = genes,
  #cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
