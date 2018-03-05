library(tximport)

args = commandArgs(trailingOnly=TRUE)

if (length(args)>0){
  print('Using input arguments')
  sams = args[1:length(args)-1]
  out_dir = args[length(args)]
} else {
  print('Using test arguments')
  sams = c('/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/ILS38024-PT1-DS1_S1/salmon/quant.sf', '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/Karpas299_S4/salmon/quant.sf', '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/RT4_S2/salmon/quant.sf', '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/RT112_S3/salmon/quant.sf')
  out_dir = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/2017-09-21_bcbio_rnaseq/expression'
}

extract_sam_name = function(path) {
  path_split = strsplit(path, split = '/')
  path_split = path_split[[1]]
  
  for ( i in 1:length(path_split)) {
    if (path_split[i] == "final"){
      return (path_split[i+1])
    }
  }
}

tx2gene = read.csv('/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/2017-09-21_bcbio_rnaseq/tx2gene.csv')

gene_counts=tximport(files = sams, type = "salmon", tx2gene = tx2gene, txOut=FALSE)
gene_counts = round(gene_counts$counts)

gene_tpm=tximport(files = sams, type = "salmon", tx2gene = tx2gene, txOut=FALSE, countsFromAbundance = "lengthScaledTPM")
gene_tpm = round(gene_tpm$counts)

isoform_tpm=tximport(files = sams, type = "salmon", tx2gene = tx2gene, txOut=TRUE, countsFromAbundance = "lengthScaledTPM")
isoform_tpm = round(isoform_tpm$counts)

names = c()
for (sam in sams){
  
  names = c(names, extract_sam_name(sam))
}
  
colnames(gene_counts)=names
colnames(gene_tpm)=names
colnames(isoform_tpm)=names

write.csv(gene_counts, file = file.path(out_dir, 'gene_counts.csv'))
write.csv(gene_tpm, file = file.path(out_dir, 'gene_tpm.csv'))
write.csv(isoform_tpm, file = file.path(out_dir, 'isoform_tpm.csv'))