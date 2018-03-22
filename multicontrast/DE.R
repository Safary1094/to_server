library(pheatmap)
library(scales)
library(gplots)
library(bcbioRNASeq)
print('bcbioRNASeq loaded')

args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  project_path = args[1]
  outputDir = args[2]
  c1 = args[3]
  c2 = args[4]
} else {
  print('Warning! No input parameters. Using defaults')
  project_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/'
  outputDir = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis/g1_g2'
  c1 = 'g1'
  c2 = 'g2'
}

#real_project_path = '/ngs/oncology/Analysis/dev/Dev_0479_HiSeq2500_STAT3_Targeted_RNASeq_set1_pool4/bcbio_rnaseq/final'

print('Reading project')
print(project_path)
bcb = loadRNASeq(uploadDir = project_path, interestingGroups = c("group"))
dds <- bcbio(bcb, "DESeqDataSet")

design(dds) = ~group
dds <- DESeq(dds)
rld <- rlog(dds)

#contrast =c("group", "DMSO_8h",  "Rosiglitazone_8h" )
contrast =c("group", c1,  c2 )
alpha = 0.01

resUnshrunken = results(dds, contrast = contrast, alpha = alpha)
resShrunken = lfcShrink(dds = dds, contrast = contrast, res = resUnshrunken)

res <- resShrunken

# prepare data
p =-log10(res@listData[["padj"]])
lfc = res@listData[["log2FoldChange"]]
baseMean = round(res@listData[["baseMean"]])
gene_id = res@rownames
# HUGO = raw_dataset$HUGO
# is_key = raw_dataset$is_key

# make dataframe
data = data.frame(gene_id, p, lfc, baseMean)
#, HUGO, is_key)


print('Writing DE result')
print(file.path(outputDir, 'RNA_DE.csv'))
write.csv(data, file.path(outputDir, 'RNA_DE.csv'), row.names=FALSE)