library(pheatmap)
library(scales)
library(gplots)
library(DESeq2)
print('DEseq2 loaded')
library(bcbioRNASeq)
print('bcbioRNASeq loaded')

args = commandArgs(trailingOnly=TRUE)

real_project_path = args[1]
#real_project_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/'

print('Reading project')
print(real_project_path)

interest_gr = c("sampleName")
bcb = loadRNASeq(uploadDir = real_project_path, interestingGroups = interest_gr)

raw_dataset = paste(bcb@metadata$projectDir, 'combined.counts', sep='/')
raw_dataset = read.csv(file=raw_dataset, header=TRUE, sep="\t")

dds <- bcbio(bcb, "DESeqDataSet")

design(dds) = ~sizeFactor+group
dds <- DESeq(dds)
rld <- rlog(dds)

contrast =c("group", "g1",  "g2" )
alpha = 0.01

resUnshrunken = results(dds, contrast = contrast, alpha = alpha)

resShrunken = lfcShrink(dds = dds, contrast = contrast, res = resUnshrunken)

res <- resShrunken

# prepare data
p =-log10(res@listData[["padj"]])
lfc = res@listData[["log2FoldChange"]]
baseMean = round(res@listData[["baseMean"]])
gene_names = res@rownames
HUGO = raw_dataset$HUGO
is_key = raw_dataset$is_key
# make dataframe
data = data.frame(p, lfc, baseMean, gene_names, HUGO, is_key)

DEout_path = args[2]
#DEout_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNA_DE.csv'

print('Writing DE result')
print(DEout_path)
write.csv(data, file = DEout_path)

# HeatMap

results <- as.data.frame(res)
results = cbind(results, data$HUGO)
results = cbind(results, data$is_key)

counts <- assay(rld)
alpha <- metadata(res)[["alpha"]]
lfc = 0
title = TRUE

results <- results %>%
  as.data.frame() %>%
  camel(strict = FALSE) %>%
  # Keep genes that pass alpha cutoff
  .[!is.na(.[["padj"]]), , drop = FALSE] %>%
  .[.[["padj"]] < alpha, , drop = FALSE] %>%
  # Keep genes that pass log2 fold change cutoff
  .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
  .[.[["log2FoldChange"]] > lfc | .[["log2FoldChange"]] < -lfc, , drop = FALSE] %>%
  .[.[["dataIsKey"]]== "True", , drop = FALSE]

genes <- rownames(results)

if (!is.null(genes)) {
  if (!all(genes %in% rownames(counts))) {
    stop(paste(
      "Genes missing from counts matrix:",
      toString(setdiff(genes, rownames(counts)))),
      call. = FALSE)
  }
  counts <- counts %>%
    .[rownames(.) %in% genes, , drop = FALSE]
} else {
  # Remove zero counts
  counts <- counts %>%
    .[rowSums(.) > 0, , drop = FALSE]
}

################################################################333
scale_rows = function (x) 
{
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}
scale_mat = function (mat, scale) 
{
  if (!(scale %in% c("none", "row", "column"))) {
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

counts_sc = scale_mat(counts, "row")

distfunc <- function(x) dist(x,method="euclidean")
hclustfunc <- function(x) hclust(x, method="ward.D")
hh<-heatmap.2(as.matrix(counts_sc),dendrogram="row",
              trace="none", margin=c(8,9), 
              hclust=hclustfunc,distfun=distfunc)

data_sorted <- counts_sc[match(rev(labels(hh$rowDendrogram)), rownames(counts_sc)), ]
aux_sorted <- results[match(rev(labels(hh$rowDendrogram)), rownames(counts_sc)), ]

data_heatmap = cbind(as.data.frame(data_sorted), aux_sorted$dataHUGO)

HMout_path = args[3]
#HMout_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNA_HM.csv'

print('Writing HM result')
print(HMout_path)
write.csv(data_heatmap, HMout_path)

