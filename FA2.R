#library("RDAVIDWebService")

# library("clusterProfiler")
# library("org.Hs.eg.db")
# library("topGO")
library("pathview")

print("FA libraries loaded")

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 2) {
  pw_path = args[1]
  genes_obj_path = args[2]
  proj_folder = args[3]
} else {
  print('Warning! No input parameters. Using defaults')
  pw_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis/RNA_PW.csv'
  genes_obj_path = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis/genes_obj.csv'
  proj_folder='/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/work/postproc/RNAanalysis'
}

pw = read.csv(pw_path)
genes_data = read.csv(genes_obj_path)

genes_obj = genes_data$V1

pathways <- pw$ID
names(genes_obj) = as.character(genes_data$X)
genes_obj = sort(genes_obj, decreasing = TRUE)

for (pw_id in pathways)
{
  a=pathview(gene.data  = genes_obj,
             pathway.id = pw_id,
             species    = "hsa",
             limit = list(gene = ceiling(max(abs(genes_obj))), cpd = 1), kegg.dir = proj_folder)
  
  gene=a$plot.data.gene
  comp=a$plot.data.cpd
  
  out_path = paste(proj_folder, '/', pw_id,'_pathway.csv', sep='')
  print(out_path)
  write.csv(rbind(gene, comp), out_path)
  
}




