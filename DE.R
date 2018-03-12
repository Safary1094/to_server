library(pheatmap)
library(scales)
library(gplots)
library(bcbioRNASeq)
print('bcbioRNASeq loaded')

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 2) {
  real_project_path = args[1]
  outputDir = args[2]
} else {
  print('Warning! No input parameters. Using defaults')
  real_project_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/'
  outputDir = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis'
}

#real_project_path = '/ngs/oncology/Analysis/dev/Dev_0479_HiSeq2500_STAT3_Targeted_RNASeq_set1_pool4/bcbio_rnaseq/final'

print('Reading project')
print(real_project_path)

iGroups <- c("group")

bcb = loadRNASeq(uploadDir = real_project_path, interestingGroups = iGroups)

###############################################################

# Directory paths (has to be inside the project directory)
#outputDir = "/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc"

if (!dir.exists(outputDir)) dir.create(outputDir)
dataDir = file.path(outputDir)
if (!dir.exists(dataDir)) dir.create(dataDir)

sampleMetadata(bcb)
rawCounts = counts(bcb, normalized = FALSE)
normalizedCounts = counts(bcb, normalized = TRUE)
rlog = counts(bcb, normalized = "rlog")
vst = counts(bcb, normalized = "vst")
tpm = tpm(bcb)

# data for Variance stabilization section + several additional files
writeCounts(rawCounts, normalizedCounts, tpm, rlog, vst, dir = dataDir, gzip=FALSE)


#######
pca=plotPCA(bcb, returnData = T)
library(DEGreport)
dds=bcb@bcbio@listData[["DESeqDataSet"]]
res = degCovariates(log2(counts(dds)+0.5),
                    colData(dds))
res2=plotPCACovariates(bcb)
corMatrix=res2$corMatrix

# data for PCA and PCA Covariates
writeCounts(corMatrix, pca, dir = dataDir, gzip=FALSE)


#######
CV=FALSE
log = "xy"
px = mcols(dds)$baseMean
sel = (px>0)
px = px[sel]
f = if (CV) sqrt else I

py = f(mcols(dds, use.names = T)$dispGeneEst[sel])
gene.names = row.names(mcols(dds, use.names = T))
gene.names = gene.names[sel]
ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)

if (!is.null(dispersions(dds))) {
  py2 = f(dispersions(dds)[sel])
}

if (!is.null(mcols(dds)$dispFit)) {
  py3 = f(mcols(dds)$dispFit[sel])
}

gene.est = as.data.frame(cbind(px, pmax(py,ymin), gene.names))
gene.fitted = as.data.frame(cbind(px, py2, gene.names))
gene.final = as.data.frame(cbind(px, py3, gene.names))

#data for the Dispersion section
writeCounts(gene.est, gene.fitted, gene.final, dir = dataDir, gzip=FALSE)

####################################################

dds <- bcbio(bcb, "DESeqDataSet")

design(dds) = ~group
dds <- DESeq(dds)
rld <- rlog(dds)

#contrast =c("group", "DMSO_8h",  "Rosiglitazone_8h" )
contrast =c("group", "g1",  "g-1" )
alpha = 0.01

resUnshrunken = results(dds, contrast = contrast, alpha = alpha)
resShrunken = lfcShrink(dds = dds, contrast = contrast, res = resUnshrunken)

res <- resShrunken

# prepare data
p =-log10(res@listData[["padj"]])
lfc = res@listData[["log2FoldChange"]]
baseMean = round(res@listData[["baseMean"]])
gene_names = res@rownames
# HUGO = raw_dataset$HUGO
# is_key = raw_dataset$is_key

# make dataframe
data = data.frame(p, lfc, baseMean, gene_names)
#, HUGO, is_key)


print('Writing DE result')
print(file.path(outputDir, 'RNA_DE.csv'))
write.csv(data, file.path(outputDir, 'RNA_DE.csv'))