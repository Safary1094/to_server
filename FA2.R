#library("RDAVIDWebService")

# library("clusterProfiler")
# library("org.Hs.eg.db")
# library("topGO")
library("pathview")

print("FA libraries loaded")

args = commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  dir_path = args[1]
} else {
  print('Warning! No input parameters. Using defaults')
  dir_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis/group_g1_vs_g2/'
}

genes_obj_path = file.path(dir_path, 'genes_obj.csv')
pathway_table_path = file.path(dir_path, 'pathway_table.csv')

pathway_table = read.csv(pathway_table_path)
genes_data = read.csv(genes_obj_path)

genes_obj = genes_data$V1

pathways <- pathway_table$ID
names(genes_obj) = as.character(genes_data$X)
genes_obj = sort(genes_obj, decreasing = TRUE)

for (pw_id in pathways)
{
  a=pathview(gene.data  = genes_obj,
             pathway.id = pw_id,
             species    = "hsa",
             limit = list(gene = ceiling(max(abs(genes_obj))), cpd = 1), kegg.dir = dir_path)
  
  gene=a$plot.data.gene
  comp=a$plot.data.cpd
  graph_file_name = paste(pw_id,'_pathway.csv', sep='')
  out_path = file.path(dir_path, graph_file_name)
  print(out_path)
  write.csv(rbind(gene, comp), out_path)
  
}




