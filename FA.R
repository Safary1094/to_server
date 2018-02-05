#library("RDAVIDWebService")

library("clusterProfiler")
library("org.Hs.eg.db")
library("topGO")
library("pathview")

print("FA libraries loaded")

args = commandArgs(trailingOnly=TRUE)



#input_path = args[1]

input_path = '~/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNA_DE.csv'

raw_dataset = read.csv(file=input_path, header=TRUE)

key_data = raw_dataset[raw_dataset$is_key == 'True',]

genes_obj = key_data[, 3]
names(genes_obj) = as.character(key_data[,5])
genes_obj = sort(genes_obj, decreasing = TRUE)

converted_names = bitr(as.character(key_data[,5]), fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
converted_names = converted_names[!duplicated(converted_names$ENSEMBL), ]
filter_key_data = key_data[is.element(key_data$gene_names,converted_names$ENSEMBL),]

filter_key_data["ENTREZid"] = converted_names$ENTREZID

genes_obj = filter_key_data[, 3]
names(genes_obj) = as.character(filter_key_data$ENTREZid)
genes_obj = sort(genes_obj, decreasing = TRUE)

gseaKEGG <- gseKEGG(
  geneList = genes_obj,
  organism = tolower("hsa"),
  nPerm = 1000,
  minGSSize = 40,
  pvalueCutoff = 0.35,
  verbose = FALSE)

gseaKEGGSummary <- slot(gseaKEGG, "result")

pathways <- gseaKEGGSummary$ID

gene_group_ens = strsplit(gseaKEGGSummary$core_enrichment[1], "/")

bitr(as.character(gene_group_ens[[1]]), fromType = "ENTREZID", toType = "ENSEMBL", OrgDb="org.Hs.eg.db")

# write files for each founded paths
#proj_folder = args[2]
proj_folder='/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/'
for (pw_id in pathways)
{
  a=pathview(gene.data  = genes_obj,
             pathway.id = pw_id,
             species    = "hsa",
             limit = list(gene = ceiling(max(abs(genes_obj))), cpd = 1))
  
  
  file.rename(paste(pw_id, '.xml', sep=''), paste(proj_folder, pw_id, '.xml', sep=''))
  
  gene=a$plot.data.gene
  comp=a$plot.data.cpd
  
  out_path = paste(proj_folder, pw_id,'_pathway.csv', sep='')
  print(out_path)
  write.csv(rbind(gene, comp), out_path)
  
}

write.csv(gseaKEGGSummary, paste(proj_folder, '/RNA_PW.csv', sep=''))


