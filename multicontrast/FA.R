#library("RDAVIDWebService")

library("clusterProfiler")
library("org.Hs.eg.db")
library("topGO")
library("pathview")

print("FA libraries loaded")

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 2) {
  input_path = args[1]
  proj_folder = args[2]
} else {
  print('Warning! No input parameters. Using defaults')
  input_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis/de_gene_key.csv'
  proj_folder='/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis'
}

data = read.csv(input_path)

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
  minGSSize = 40,
  pvalueCutoff = 0.50,
  verbose = FALSE)

gseaKEGGSummary <- slot(gseaKEGG, "result")

pathways <- gseaKEGGSummary$ID

gene_group_ens = strsplit(gseaKEGGSummary$core_enrichment[1], "/")

bitr(as.character(gene_group_ens[[1]]), fromType = "ENTREZID", toType = "ENSEMBL", OrgDb="org.Hs.eg.db")

write.csv(gseaKEGGSummary, paste(proj_folder, '/RNA_PW.csv', sep=''))
write.csv(as.matrix(genes_obj), paste(proj_folder, '/genes_obj.csv', sep=''))
# write files for each founded paths



for (pw_id in pathways)
{
  a=pathview(gene.data  = genes_obj,
             pathway.id = pw_id,
             species    = "hsa",
             limit = list(gene = ceiling(max(abs(genes_obj))), cpd = 1), kegg.dir = proj_folder)
  
  #file.rename(paste(pw_id, '.xml', sep=''), paste(proj_folder, '/', pw_id, '.xml', sep=''))
  
  gene=a$plot.data.gene
  comp=a$plot.data.cpd
  
  out_path = paste(proj_folder, '/', pw_id,'_pathway.csv', sep='')
  print(out_path)
  write.csv(rbind(gene, comp), out_path)
  
}




