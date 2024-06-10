library("dplyr")
library("stringr")

chromosome_length <- read.table("/genome/mm10/mm10.chrom.sizes", header=FALSE)
chr_length <- chromosome_length$V2
names(chr_length) <- chromosome_length$V1

gene_tss <- read.table("/annot/gencode.vM23.gene.tssUpDn1k.bed")

gene_tss %>% mutate(V9=V2-500000, V10=V3+500000) %>% mutate(V11=case_when(V9<0 ~ 0, .default=V9), V12=case_when(V10>chr_length[V1] ~ chr_length[V1], .default=V10)) -> gene_tss

gene_tss %>% select(V1, V11, V12, V4, V5, V6, V7, V8) -> gene_tss

write.table(gene_tss, file="/annot/gencode.vM23.gene.tssUpDn500k.bed", quote=FALSE, row.names=FALSE, col.name=FALSE, sep="\t")
