library(DESeq2)

args <- commandArgs(TRUE)
csvFilePath = args[1]
#csvFilePath = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/pca_data.csv'
pcaOutputFpath = args[2]
#geneCombinedCounts = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/2017-09-21_bcbio_rnaseq/expression/gene_counts.csv'
geneCombinedCounts = args[3]

meta <- read.csv(csvFilePath, stringsAsFactors=F, header=T)
rownames(meta) <- meta$description
if(!is.element("condition", colnames(meta))){{
  meta <- data.frame(meta, condition1=rep("group1", nrow(meta)), condition2=meta$description)
}} else {{
  meta <- data.frame(meta, condition1=meta$condition, condition2=meta$condition)
}}

metax <- data.frame(name=as.character(meta$description), condition=as.character(meta$condition2))
rownames(metax) <-  metax$name

countsx <- read.csv(geneCombinedCounts, header=F, stringsAsFactors = F)
rownames(countsx) <- countsx[,1]
colnames(countsx) <- countsx[1,]
countsx <- countsx[-1,-1]
countsx <- countsx[, rownames(metax)]
counts <- apply(countsx, 2, as.numeric)
rownames(counts) <- rownames(countsx)

dds <- DESeqDataSetFromMatrix(countData = counts, colData=metax, design= ~ condition)
rld <- rlog(dds)

data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE, ntop=1000)
percentVar <- round(100 * attr(data, "percentVar"))
write.table(data, file=pcaOutputFpath, sep="\t")
line = paste0("#Variance:", percentVar[1], ",", percentVar[2])
write(line, file=pcaOutputFpath, append=TRUE)