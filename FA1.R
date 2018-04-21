#library("RDAVIDWebService")

library("clusterProfiler")
library("org.Hs.eg.db")
library("topGO")
#library("pathview")

print("FA libraries loaded")

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 2) {
  dir_path = args[1]
  de_genes = args[2]
} else {
  print('Warning! No input parameters. Using defaults')
  de_genes = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis/group_g3_vs_g4/de_gene_all.csv'
  dir_path='/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis/group_g3_vs_g4/'
}

data = read.csv(de_genes)

converted_names = bitr(as.character(data$gene_id), fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
converted_names = converted_names[!duplicated(converted_names$ENSEMBL), ]
data = data[is.element(data$gene_id, converted_names$ENSEMBL),]

data["ENTREZid"] = converted_names$ENTREZID

genes_obj = data$p
names(genes_obj) = as.character(data$ENTREZid)
genes_obj = sort(genes_obj, decreasing = TRUE)

gseaKEGG <- gseKEGG(
  geneList = genes_obj,
  organism = tolower("hsa"),
  nPerm = 1000,
  minGSSize = 100,
  pvalueCutoff = 0.05,
  verbose = FALSE)

gseaKEGGSummary <- slot(gseaKEGG, "result")
print(gseaKEGGSummary$ID)

pathways <- gseaKEGGSummary$ID

gene_group_ens = strsplit(gseaKEGGSummary$core_enrichment[1], "/")

bitr(as.character(gene_group_ens[[1]]), fromType = "ENTREZID", toType = "ENSEMBL", OrgDb="org.Hs.eg.db")

write.csv(gseaKEGGSummary, file.path(dir_path, 'pathway_table.csv'), row.names=FALSE)
write.csv(as.matrix(genes_obj), file.path(dir_path, 'genes_obj.csv'))