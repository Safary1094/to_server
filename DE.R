library(pheatmap)
library(scales)
library(gplots)
library(bcbioRNASeq)
print('bcbioRNASeq loaded')

args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  project_path = args[1]
  outputDir = args[2]
  contrast_path = args[3]
} else {
  print('Warning! No input parameters. Using defaults')
  project_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/'
  outputDir = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis/'
  contrast_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/2017-09-21_bcbio_rnaseq/contrasts.csv'
}

#real_project_path = '/ngs/oncology/Analysis/dev/Dev_0479_HiSeq2500_STAT3_Targeted_RNASeq_set1_pool4/bcbio_rnaseq/final'

print('Reading project')
print(project_path)
bcb = loadRNASeq(uploadDir = project_path, interestingGroups = c("group"))
dds <- bcbio(bcb, "DESeqDataSet")

design(dds) = ~group
dds <- DESeq(dds)

contrast =list()
all_contrasts = as.matrix(read.csv(contrast_path, header = FALSE))
for (i in 1: nrow(all_contrasts))
{
  cont = unname(all_contrasts[i,])
  contrast[[i]] = c("group", cont[1], cont[2] )
}
 
comp_res=degComps(dds, contrast = contrast)

i=1
for (res in comp_res)
{
  d = res$shrunken
  # prepare data
  p =-log10(d$padj)
  lfc = d$log2FoldChange
  baseMean = round(d$baseMean)
  gene_id = row.names(d)
  
  # make dataframe
  data = data.frame(gene_id, p, lfc, baseMean)
  
  contrast_name = names(comp_res[i])
  
  writetDir = file.path(outputDir, contrast_name)
  if (!dir.exists(writetDir)) dir.create(writetDir)
  print('Writing DE result')
  print(file.path(writetDir,'RNA_DE.csv'))
  write.csv(data, file.path(writetDir, 'RNA_DE.csv'), row.names=FALSE)

  i = i+1
}