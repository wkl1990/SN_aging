suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
library("biomaRt")

# Load RNA data to get all gene names
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
all_genes <- rownames(RNA)

# Basic function to convert mouse to human gene names
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
mouse2human <- getLDS(attributes = c("ensembl_gene_id","external_gene_name", "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_orthology_type"),
       filters = "external_gene_name", values = all_genes,
       mart = mouse,
       attributesL = c("ensembl_gene_id","external_gene_name", "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type"), 
       martL = human)
write.csv(mouse2human, file="./rds/ATAC/after_integra/GRN/SCENT/VShuman/mouse2human_gene_conversion.csv", row.names=FALSE)


# reciprocal cCREs
all_peak <- read.table("/projects/ps-renlab2/kaw033/WK_aging/peak_calling/after_integra/final/SN_integra.final.peak.srt.bed")
colnames(all_peak) <- c("chrom","start","end","peak")

reciprocal_peak <- read.table("/projects/ps-renlab2/kaw033/WK_aging/ldsc/after_integra/celltype/raw/SN_integra.final.peak.srt.reciprocalToHg38.bed")
colnames(reciprocal_peak) <- c("chrom","start","end","Mpeak")
length(intersect(all_peak$peak, reciprocal_peak$Mpeak))
reciprocal_peak %>% mutate(Hpeak = paste0(chrom, ":", start, "-", end)) -> reciprocal_peak

# get human pdc from Yang
human_peak <- read.table("rds/ATAC/after_integra/GRN/SCENT/VShuman/YangLi_cCRE/Table S6 – List of cCREs in bed format")
colnames(human_peak) <- c("chrom","start","end","cCRE")
human_peak$peak <- gsub("cCRE", "peak", human_peak$cCRE)

human_pdc <- read.table("rds/ATAC/after_integra/GRN/SCENT/VShuman/YangLi_cCRE/Table S12 - Summary of gene-cCRE correlations.txt", header=TRUE)
human_pdc %>% left_join(human_peak, by=c("distal_cCRE"="peak")) %>% dplyr::select(chrom, start, end, distal_cCRE, gene, proximal_cCRE) -> human_pdc_bed

# get mouse pdc
ATAC_subclass_forCicero <- read.csv("rds/ATAC/after_integra/GRN/cicero/subclass_list.csv", header=FALSE)
scent_pdc_all <- read.table(file="/projects/ps-renlab2/kaw033/WK_aging/rds/ATAC/after_integra/GRN/SCENT/peak_gene_link/SN_all.SCENT.bedpe")
colnames(scent_pdc_all) <- c("chr1","start1","end1","chrom","start","end","pdc","score","V9","V10")

scent_pdc_all %>% tidyr::separate(pdc, into=c("gene","peak"), sep="[|]", remove=FALSE) -> scent_pdc_all
scent_pdc_all$Hpeak <- reciprocal_peak$Hpeak[match(scent_pdc_all$peak, reciprocal_peak$Mpeak)]
mouse2human %>% dplyr::select(Gene.name, Human.gene.name) %>% rename(gene=Gene.name, Hgene=Human.gene.name) -> mouse2human_rename
scent_pdc_all_gene <- scent_pdc_all %>% left_join(mouse2human_rename, by=join_by("gene"=="gene"))
#scent_pdc_all$Hgene <- mouse2human$Human.gene.name[match(scent_pdc_all$gene, mouse2human$Gene.name)]
scent_pdc_all_gene$Hgene2UP <- toupper(scent_pdc_all_gene$gene)
scent_pdc_all_gene %>% mutate(Hgene2UPvsHgene = ifelse(Hgene2UP == Hgene, TRUE, FALSE), Hgene_final=ifelse(is.na(Hgene) | Hgene=="", Hgene2UP, Hgene)) -> scent_pdc_all_gene
scent_pdc_all_gene %>% filter(Hgene2UPvsHgene!=TRUE) %>% dplyr::select(gene, Hgene, Hgene2UP, Hgene_final) %>% distinct() -> scent_pdc_gene_conversion_check

scent_pdc_all_gene %>% filter(is.na(Hpeak)) -> scent_pdc_noReciprocal_peak
scent_pdc_all_gene %>% filter(!is.na(Hpeak)) -> scent_pdc_withReciprocal_peak

scent_pdc_withReciprocal_peak %>% dplyr::select(Hpeak, Hgene_final, Hgene2UP, gene, pdc) %>% tidyr::separate(Hpeak, into=c("chrom", "start", "end"), sep="[:-]", remove=FALSE) %>% dplyr::select(chrom, start, end, Hpeak, Hgene_final,Hgene2UP, gene, pdc) %>% distinct() -> scent_pdc_withReciprocal_peak_bed

bedtoolsr::bt.intersect(scent_pdc_withReciprocal_peak_bed, human_pdc_bed, wao=TRUE) %>% filter(V9 != ".") %>% dplyr::select(V5, V6, V7, V13) -> scent_human_overlap_gene_check
write.csv(scent_human_overlap_gene_check, file="./rds/ATAC/after_integra/GRN/SCENT/VShuman/scent_human_overlap_gene_check.csv", row.names=FALSE)

bedtoolsr::bt.intersect(scent_pdc_withReciprocal_peak_bed, human_pdc_bed, wao=TRUE) %>% filter(V9 != ".") %>% mutate(match=ifelse(V5 == V13 | V6 == V13, TRUE, FALSE)) %>% filter(match==TRUE) %>% 
       mutate(human_peak_coord=paste0(V9, ":", V10, "-", V11)) %>% dplyr::select(V8, V4, V12, human_peak_coord, V13) %>% rename(PDC=V8, peak_in_hg38=V4, human_ccre=V12, human_gene=V13) %>% distinct() -> scent_human_overlap_final
write.csv(scent_human_overlap_final, file="./rds/ATAC/after_integra/GRN/SCENT/VShuman/scent_human_overlap_final.csv", row.names=FALSE)


# for each cell type
scent_pdc_list <- list()
scent_pdc_overlap_list <- list()
for (subclass_label_id in ATAC_subclass_forCicero$V1) {
  subclass_name <- gsub(" ", "_", subclass_label_id)
  file <- file.path("/projects/ps-renlab2/kaw033/WK_aging/", "rds/ATAC/after_integra/GRN/SCENT/peak_gene_link/", paste0(subclass_name, ".SCENT.bedpe"))
  scent_pdc_celltype<- read.table(file)
  colnames(scent_pdc_celltype) <- c("chr1","start1","end1","chrom","start","end","pdc","score","V9","V10")
  scent_pdc_celltype %>% tidyr::separate(pdc, into=c("gene","peak"), sep="[|]", remove=FALSE) -> scent_pdc_celltype

  scent_pdc_celltype$Hpeak <- reciprocal_peak$Hpeak[match(scent_pdc_celltype$peak, reciprocal_peak$Mpeak)]

  scent_pdc_celltype_gene <- scent_pdc_celltype %>% left_join(mouse2human_rename, by=join_by("gene"=="gene"))
  #scent_pdc_celltype$Hgene <- mouse2human$Human.gene.name[match(scent_pdc_celltype$gene, mouse2human$Gene.name)]
  scent_pdc_celltype_gene$Hgene2UP <- toupper(scent_pdc_celltype_gene$gene)
  scent_pdc_celltype_gene %>% mutate(Hgene2UPvsHgene = ifelse(Hgene2UP == Hgene, TRUE, FALSE), Hgene_final=ifelse(is.na(Hgene) | Hgene=="", Hgene2UP, Hgene)) -> scent_pdc_celltype_gene
  scent_pdc_list[[subclass_name]] <- scent_pdc_celltype

  scent_pdc_celltype_gene %>% filter(is.na(Hpeak)) -> scent_pdc_celltype_noReciprocal_peak
  scent_pdc_celltype_gene %>% filter(!is.na(Hpeak)) -> scent_pdc_celltype_withReciprocal_peak

  scent_pdc_celltype_withReciprocal_peak %>% dplyr::select(Hpeak, Hgene_final, Hgene2UP, gene, pdc) %>% tidyr::separate(Hpeak, into=c("chrom", "start", "end"), sep="[:-]", remove=FALSE) %>% dplyr::select(chrom, start, end, Hpeak, Hgene_final,Hgene2UP, gene, pdc) %>% distinct() -> scent_pdc_celltype_withReciprocal_peak_bed

  bedtoolsr::bt.intersect(scent_pdc_celltype_withReciprocal_peak_bed, human_pdc_bed, wao=TRUE) %>% filter(V9 != ".") %>% mutate(match=ifelse(V5 == V13 | V6 == V13, TRUE, FALSE)) %>% filter(match==TRUE) %>% 
       mutate(human_peak_coord=paste0(V9, ":", V10, "-", V11)) %>% dplyr::select(V8, V4, V12, human_peak_coord, V13) %>% rename(PDC=V8, peak_in_hg38=V4, human_ccre=V12, human_gene=V13) %>% distinct() -> scent_celltype_human_overlap_final
  write.csv(scent_celltype_human_overlap_final, file=paste0("./rds/ATAC/after_integra/GRN/SCENT/VShuman/scent_human_overlap_final_", subclass_name, ".csv"), row.names=FALSE)
  scent_pdc_overlap_list[[subclass_name]] <- scent_celltype_human_overlap_final
}
