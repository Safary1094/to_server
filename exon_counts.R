library(Rsubread)

args = commandArgs(trailingOnly=TRUE)

if (length(args)>0){
  print('Using input arguments')
  bams = args[1:(length(args)-2)]
  out_dir = args[length(args)-1]
  ann = args[length(args)]
} else {
  print('Using test arguments')
  bams = c('/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/ILS38024-PT1-DS1_S1/ILS38024-PT1-DS1_S1-ready.bam', '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/Karpas299_S4/Karpas299_S4-ready.bam', '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/RT4_S2/RT4_S2-ready.bam', '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/RT112_S3/RT112_S3-ready.bam')
  out_dir = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/Dev_0406/final/2017-09-21_bcbio_rnaseq/expression'
  ann = '/home/alexey/Downloads/ref-transcripts.dexseq.gff3'
}

extract_bam_name = function(path) {
  path_split = strsplit(path, split = '/')
  path_split = path_split[[1]]
  
  for ( i in 1:length(path_split)) {
    if (path_split[i] == "final"){
      return (path_split[i+1])
    }
  }
}

print('bams:')
print(bams)
print('ann:')
print(ann)

res_exons <- featureCounts(files = bams, annot.ext = ann, isGTFAnnotationFile = T,  useMetaFeatures = F, allowMultiOverlap = T, nthreads = 4, countMultiMappingReads = F, GTF.featureType = 'exonic_part', fraction = 0.9, isPairedEnd = T)

data = res_exons$counts

names = c()
for (bam in bams){
  names = c(names, extract_bam_name(bam))
}

colnames(data)=names

print('out_dir:')
print(out_dir)
write.csv(data, file.path(out_dir, 'exon_counts.csv'))


 