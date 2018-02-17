library(Rsubread)

ann = '/ngs/reference_data/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.gtf'
#ann = '/home/alexey/Downloads/ref-transcripts.cancer.trx.gtf'

bam = '/ngs/oncology/Analysis/external/EXT_087_Novogene_CDK9i_RNASeq/rnaseq/final/TR4h_Repeat1'
#bam = '/home/alexey/ngs/NGS_Reporting_TestData/data/bcbio_postproc/bams/syn3-normal-ready.bam'

outp = '/ngs/usr/safary/exon.csv'

#GTF.attrType = "exon_id",
res_exons <- featureCounts(files = bam, annot.ext = ann, isGTFAnnotationFile = T,  useMetaFeatures = F)

data = cbind(res_exons$counts, res_exons$annotation$Start, res_exons$annotation$End)

n = as.data.frame(res_exons$annotation$GeneID)
c = as.data.frame(res_exons$counts)
s = as.data.frame(res_exons$annotation$Start)
e = as.data.frame(res_exons$annotation$End)

a = data.frame(n,c,s,e)

colnames(a)[2] <- "counts"

write.csv(a, outp)
