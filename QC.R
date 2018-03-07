###__Quality_control__###

library(bcbioRNASeq)

# Load bcbioRNASeq object
args = commandArgs(trailingOnly=TRUE)

project_path = args[1]
#project_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/'

bcb = loadRNASeq(project_path)

# Directory paths (has to be inside the project directory)
outputDir = args[2]
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