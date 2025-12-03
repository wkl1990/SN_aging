suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("monocle3"))

# check PD KEGG genes
# PD kegg pathway
Autochr.genetable_enzid <- read.csv(file="rds/RNA/preprocess/Autochr.genetable.entrezid.csv", row.names=1)
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)
RNA_subclass_forDEG <- intersect(rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell)>=100], rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0])


PD_KEGG_monocle3 <- list()
for (subclass_label_id in RNA_subclass_forDEG) {
  subclass_name <- gsub(" ", "_", subclass_label_id)
  GO_KEGG_file <- paste0("rds/RNA/diff/after_integra/monocle3/subclass/", subclass_name, ".deg.GO_KEGG.csv")
  if (file.exists(GO_KEGG_file)) {
    GO_KEGG <- read.csv(paste0("rds/RNA/diff/after_integra/monocle3/subclass/", subclass_name, ".deg.GO_KEGG.csv"), row.names=1)
    GO_KEGG_PD <- GO_KEGG[grepl("Parkinson", GO_KEGG$Description), ]
    if(nrow(GO_KEGG_PD) != 0) {
      PD_KEGG_monocle3[[subclass_label_id]] <- GO_KEGG_PD
    }
  }
}

PD_KEGG_monocle3_tbl <- do.call(rbind,PD_KEGG_monocle3) 
PD_KEGG_monocle3_tbl %>% select(group, p.adjust, geneID) %>% filter(group=="up") %>% tibble::rownames_to_column("Description") %>% mutate(Description=gsub("\\.up\\.mmu05012", "", Description)) -> PD_KEGG_monocle3_data

PD_KEGG_monocle3_genelist <- list()
for (subclass in PD_KEGG_monocle3_data$Description) {
	geneID <- PD_KEGG_monocle3_data$geneID[which(PD_KEGG_monocle3_data$Description==subclass)]
	geneIDs <- unlist(strsplit(geneID, "/"))
	gene_symbols <- Autochr.genetable_enzid$SYMBOL[match(geneIDs, Autochr.genetable_enzid$ENTREZID)]
	PD_KEGG_monocle3_genelist[[subclass]] <- gene_symbols
}

#maxlen <- max(sapply(PD_KEGG_monocle3_genelist,length))
#lapply(seq(maxlen),function(i) Reduce(intersect,lapply(PD_KEGG_monocle3_genelist,"[[",i)))
intersect(intersect(intersect(intersect(intersect(intersect(PD_KEGG_monocle3_genelist[[1]],PD_KEGG_monocle3_genelist[[2]]),PD_KEGG_monocle3_genelist[[3]]),PD_KEGG_monocle3_genelist[[4]]),PD_KEGG_monocle3_genelist[[5]]),PD_KEGG_monocle3_genelist[[6]]),PD_KEGG_monocle3_genelist[[7]])
intersect(intersect(PD_KEGG_monocle3_genelist[["166 MRN Pou3f1 C1ql4 Glut"]],PD_KEGG_monocle3_genelist[["197 SNr Six3 Gaba"]]),PD_KEGG_monocle3_genelist[["327 Oligo NN"]])

PD_KEGG_NOISeq <- list()
for (subclass_label_id in RNA_subclass_forDEG) {
  subclass_name <- gsub(" ", "_", subclass_label_id)
  GO_KEGG_file <- paste0("rds/RNA/diff/after_integra/NOISeq/subclass/", subclass_name, ".deg.GO_KEGG.csv")
  if (file.exists(GO_KEGG_file)) {
    GO_KEGG <- read.csv(paste0("rds/RNA/diff/after_integra/NOISeq/subclass/", subclass_name, ".deg.GO_KEGG.csv"), row.names=1)
    GO_KEGG_PD <- GO_KEGG[grepl("Parkinson", GO_KEGG$Description), ]
    if(nrow(GO_KEGG_PD) != 0) {
      PD_KEGG_NOISeq[[subclass_label_id]] <- GO_KEGG_PD
    }
  }
}

PD_KEGG_NOISeq_tbl <- do.call(rbind,PD_KEGG_NOISeq) 
PD_KEGG_NOISeq_tbl %>% select(group, p.adjust, geneID) %>% filter(group=="up") %>% tibble::rownames_to_column("Description") %>% mutate(Description=gsub("\\.up\\.mmu05012", "", Description)) -> PD_KEGG_NOISeq_data

PD_KEGG_NOISeq_genelist <- list()
for (subclass in PD_KEGG_NOISeq_data$Description) {
	geneID <- PD_KEGG_NOISeq_data$geneID[which(PD_KEGG_NOISeq_data$Description==subclass)]
	geneIDs <- unlist(strsplit(geneID, "/"))
	gene_symbols <- Autochr.genetable_enzid$SYMBOL[match(geneIDs, Autochr.genetable_enzid$ENTREZID)]
	PD_KEGG_NOISeq_genelist[[subclass]] <- gene_symbols
}

intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(PD_KEGG_NOISeq_genelist[[1]],PD_KEGG_NOISeq_genelist[[2]]),PD_KEGG_NOISeq_genelist[[3]]),PD_KEGG_NOISeq_genelist[[4]]),
	PD_KEGG_NOISeq_genelist[[5]]),PD_KEGG_NOISeq_genelist[[6]]),PD_KEGG_NOISeq_genelist[[7]]),PD_KEGG_NOISeq_genelist[[8]]),PD_KEGG_NOISeq_genelist[[9]]),PD_KEGG_NOISeq_genelist[[10]]),
	PD_KEGG_NOISeq_genelist[[11]]),PD_KEGG_NOISeq_genelist[[12]]),PD_KEGG_NOISeq_genelist[[13]]),PD_KEGG_NOISeq_genelist[[14]]),PD_KEGG_NOISeq_genelist[[15]]),PD_KEGG_NOISeq_genelist[[16]]),PD_KEGG_NOISeq_genelist[[17]])
intersect(intersect(intersect(intersect(intersect(PD_KEGG_NOISeq_genelist[["166 MRN Pou3f1 C1ql4 Glut"]],PD_KEGG_NOISeq_genelist[["197 SNr Six3 Gaba"]]),PD_KEGG_NOISeq_genelist[["327 Oligo NN"]]),PD_KEGG_monocle3_genelist[["166 MRN Pou3f1 C1ql4 Glut"]]),PD_KEGG_monocle3_genelist[["197 SNr Six3 Gaba"]]),PD_KEGG_monocle3_genelist[["327 Oligo NN"]])


# overlap with PD genes
expressed_gene <- function(file, express=0, cell=1) {
	express_mat <- read.csv(file, row.names=1)
	expressGene <- rownames(express_mat)[rowSums(express_mat>express)>=cell]
	return(expressGene)	
}

raregene <- read.table("human_PD/PDgenes/t_rare_gene20200820.txt", skip=1, sep="\t")
rare_genes <- unique(raregene[,1])
rare_genes <- trimws(rare_genes)
rare_genes <- sub("\\?","",rare_genes)
rare_genes <- toupper(rare_genes)
rare_genes[grepl("-", rare_genes)]
rare_genes[is.na(rare_genes)]
rare_genes[which(rare_genes=="-")]
rare_genes[which(rare_genes=="")]

expgene <- read.table("human_PD/PDgenes/t_gene_expression20200820.txt", skip=1, sep="\t", fill=TRUE)
exp_genes <- unique(expgene[,2])
exp_genes <- trimws(exp_genes)
exp_genes <- toupper(exp_genes)
exp_genes[grepl("-", exp_genes)]
exp_genes <- sub("6-SEP", "SEPT6", exp_genes)
exp_genes <- sub("5-MAR", "MARCH5", exp_genes)
exp_genes <- sub("10-SEP", "SEPT10", exp_genes)
exp_genes <- sub("4-SEP", "SEPT4", exp_genes)
exp_genes[is.na(exp_genes)]
exp_genes[which(exp_genes=="-")]
exp_genes[which(exp_genes=="")]

methgene <- read.table("human_PD/PDgenes/t_dna_methylation20200820.txt", sep="\t", skip=1)
meth_genes <- unique(methgene[,4])
meth_genes <- trimws(meth_genes)
meth_genes <- toupper(meth_genes)
meth_genes[grepl("-", meth_genes)]
meth_genes[is.na(meth_genes)]
meth_genes[which(meth_genes=="-")]
meth_genes[which(meth_genes=="")]
meth_genes <- meth_genes[which(meth_genes!="-")]
meth_genes <- meth_genes[which(meth_genes!="")]

rareVgene <- read.table("human_PD/PDgenes/t_rare_variant20200820.txt", skip=1, sep="\t", quote="", fill=TRUE)
rareV_genes <- unique(rareVgene[,10])
rareV_genes <- trimws(rareV_genes)
rareV_genes <- toupper(rareV_genes)
rareV_genes[grepl("-", rareV_genes)]
rareV_genes[grepl(";", rareV_genes)]
rareV_genes <- c(rareV_genes, "MIR12136", "FAM47E", "FAM47E-STBD1", "PINK1", "PINK1-AS")
rareV_genes <- unique(rareV_genes[!grepl(";", rareV_genes)])
rareV_genes[is.na(rareV_genes)]
rareV_genes[which(rareV_genes=="-")]
rareV_genes[which(rareV_genes=="")]

gwasgene <- read.table("human_PD/PDgenes/t_common_variant20200820.txt", skip=1, sep="\t")
gwas_genes <- unique(gwasgene[,2])
gwas_genes <- unlist(lapply(gwas_genes, function(x) str_split(x, ";|,|/")))
gwas_genes <- trimws(gwas_genes)
gwas_genes <- toupper(gwas_genes)
gwas_genes[grepl("-", gwas_genes)]
gwas_genes[grepl(";", gwas_genes)]
gwas_genes[grepl("=", gwas_genes)]
gwas_genes <- sub("=[0-9]*", "", gwas_genes)
gwas_genes[is.na(gwas_genes)]
gwas_genes[which(gwas_genes=="-")]
gwas_genes[which(gwas_genes=="")]
gwas_genes <- gwas_genes[!is.na(gwas_genes)]
gwas_genes <- gwas_genes[which(gwas_genes!="-")]
gwas_genes <- gwas_genes[which(gwas_genes!="")]

cnvgene <- read.table("human_PD/PDgenes/t_cnv20200820.txt", skip=1, sep="\t", quote="")
cnv_genes <- unique(cnvgene[,1])
cnv_genes <- trimws(cnv_genes)
cnv_genes <- toupper(cnv_genes)
cnv_genes[grepl("-", cnv_genes)]
cnv_genes[is.na(cnv_genes)]
cnv_genes[which(cnv_genes=="-")]
cnv_genes[which(cnv_genes=="")]

gene4PDgene <- read.table("human_PD/PDgenes/Gene4PD_Gene.txt", sep="\t", quote="", skip=1)
gene4PD_genes <- unique(gene4PDgene[,1])
gene4PD_genes <- trimws(gene4PD_genes)
gene4PD_genes <- toupper(gene4PD_genes)
gene4PD_genes[grepl("-", gene4PD_genes)]
gene4PD_genes[is.na(gene4PD_genes)]
gene4PD_genes[which(gene4PD_genes=="-")]
gene4PD_genes[which(gene4PD_genes=="")]

malagene <- read.table("human_PD/PDgenes/Malacards.txt", row.names=1, fill=TRUE, skip=1)
mala_genes <- unique(malagene[,1])
mala_genes <- trimws(mala_genes)
mala_genes <- toupper(mala_genes)
mala_genes[grepl("-", mala_genes)]
mala_genes[is.na(mala_genes)]
mala_genes[which(mala_genes=="-")]
mala_genes[which(mala_genes=="")]

mala_high_genes <- unique(malagene[which(malagene[,2]=="high"),1])
mala_high_genes <- trimws(mala_high_genes)
mala_high_genes <- toupper(mala_high_genes)
mala_high_genes[grepl("-", mala_high_genes)]
mala_high_genes[is.na(mala_high_genes)]
mala_high_genes[which(mala_high_genes=="-")]
mala_high_genes[which(mala_high_genes=="")]

malagene_ALS <- read.table("human_PD/PDgenes/Malacards_ALS1.txt", row.names=1, fill=TRUE, skip=1)
malaALS_genes <- unique(malagene_ALS[,1])
malaALS_genes <- trimws(malaALS_genes)
malaALS_genes <- toupper(malaALS_genes)
malaALS_genes[grepl("-", malaALS_genes)]
malaALS_genes[is.na(malaALS_genes)]
malaALS_genes[which(malaALS_genes=="-")]
malaALS_genes[which(malaALS_genes=="")]


save(rare_genes, exp_genes, meth_genes, rareV_genes, gwas_genes, cnv_genes, gene4PD_genes, mala_genes, mala_high_genes, malaALS_genes, file="./human_PD/PDgenes/PDgenes.RDS")

DEG_monocle3_PDgene_pval <- list()
DEG_monocle3_PDgene_raregenes <- list()
DEG_monocle3_PDgene_expgenes <- list()
DEG_monocle3_PDgene_methgenes <- list()
DEG_monocle3_PDgene_rareVgenes <- list()
DEG_monocle3_PDgene_gwasgenes <- list()
DEG_monocle3_PDgene_cnvgenes <- list()
DEG_monocle3_PDgene_gene4PDgenes <- list()
DEG_monocle3_PDgene_malagenes <- list()
DEG_monocle3_PDgene_malahighgenes <- list()
DEG_monocle3_PDgene_malaALSgenes <- list()

for (subclass_label_id in RNA_subclass_forDEG) {
	print(paste0(subclass_label_id, " start!"))
	# read deg
	subclass_name <- gsub(" ", "_", subclass_label_id)
	deg_file <- paste0("./rds/RNA/diff/after_integra/monocle3/subclass/", subclass_name, ".age_terms_deg.csv")
	celltype_age_terms_deg <- read.csv(deg_file, row.names=1)
	if (nrow(celltype_age_terms_deg)<1) {
		next
	}
	deg_list <- list()
	deg_list[["all"]] <- celltype_age_terms_deg %>% pull(name) %>% toupper
	deg_list[["up"]] <- celltype_age_terms_deg %>% dplyr::filter(estimate>0) %>% pull(name) %>% toupper
	deg_list[["down"]] <- celltype_age_terms_deg %>% dplyr::filter(estimate<0) %>% pull(name) %>% toupper
	express_file <- paste0("./rds/RNA/diff/after_integra/pseudo_bulk/subclass/", subclass_name, "_pseudobulk_autochr.csv")
	express_genes <- expressed_gene(express_file, express=0, cell=1)	
	pval_list <- c()
	for (type in c("all", "up", "down")) {
		for (genes in list(rare_genes, exp_genes, meth_genes, rareV_genes, gwas_genes, cnv_genes, gene4PD_genes, mala_genes, mala_high_genes, malaALS_genes)) {
			ovlp_len <- length(intersect(deg_list[[type]], genes))
			if (ovlp_len>0) {
				pval <- 1-phyper(length(intersect(deg_list[[type]], genes)), length(deg_list[[type]]), length(express_genes)-length(deg_list[[type]]), length(genes))
			} else {
				pval <- NA
			}			
			pval_list <- c(pval_list, pval)
		}
	}
	DEG_monocle3_PDgene_pval[[subclass_label_id]] <- pval_list

	ovlp_raregenes <- intersect(deg_list[["up"]], rare_genes)
	ovlp_expgenes <- intersect(deg_list[["up"]], exp_genes)
	ovlp_methgenes <- intersect(deg_list[["up"]], meth_genes)
	ovlp_rareVgenes <- intersect(deg_list[["up"]], rareV_genes)
	ovlp_gwasgenes <- intersect(deg_list[["up"]], gwas_genes)
	ovlp_cnvgenes <- intersect(deg_list[["up"]], cnv_genes)
	ovlp_gene4PDgenes <- intersect(deg_list[["up"]], gene4PD_genes)
	ovlp_malagenes <- intersect(deg_list[["up"]], mala_genes)
	ovlp_malahighgenes <- intersect(deg_list[["up"]], mala_high_genes)
	ovlp_malaALSgenes <- intersect(deg_list[["up"]], malaALS_genes)
	DEG_monocle3_PDgene_raregenes[[subclass_label_id]] <- ovlp_raregenes
	DEG_monocle3_PDgene_expgenes[[subclass_label_id]] <- ovlp_expgenes
	DEG_monocle3_PDgene_methgenes[[subclass_label_id]] <- ovlp_methgenes
	DEG_monocle3_PDgene_rareVgenes[[subclass_label_id]] <- ovlp_rareVgenes
	DEG_monocle3_PDgene_gwasgenes[[subclass_label_id]] <- ovlp_gwasgenes
	DEG_monocle3_PDgene_cnvgenes[[subclass_label_id]] <- ovlp_cnvgenes
	DEG_monocle3_PDgene_gene4PDgenes[[subclass_label_id]] <- ovlp_gene4PDgenes
	DEG_monocle3_PDgene_malagenes[[subclass_label_id]] <- ovlp_malagenes
	DEG_monocle3_PDgene_malahighgenes[[subclass_label_id]] <- ovlp_malahighgenes
	DEG_monocle3_PDgene_malaALSgenes[[subclass_label_id]] <- ovlp_malaALSgenes
}
DEG_monocle3_PDgene_pval_combined <- do.call(rbind, DEG_monocle3_PDgene_pval)
colnames(DEG_monocle3_PDgene_pval_combined) <- c("all_rare_genes", "all_exp_genes", "all_meth_genes", "all_rareV_genes", "all_gwas_genes", "all_cnv_genes", "all_gene4PD_genes", "all_mala_genes", "all_mala_high_genes", "all_malaALS_genes", 
	"up_rare_genes", "up_exp_genes", "up_meth_genes", "up_rareV_genes", "up_gwas_genes", "up_cnv_genes", "up_gene4PD_genes", "up_mala_genes", "up_mala_high_genes", "up_malaALS_genes", 
	"down_rare_genes", "down_exp_genes", "down_meth_genes", "down_rareV_genes", "down_gwas_genes", "down_cnv_genes", "down_gene4PD_genes", "down_mala_genes", "down_mala_high_genes", , "down_malaALS_genes")
write.csv(DEG_monocle3_PDgene_pval_combined, file="./human_PD/PDgenes/deg.monocle3.PDgenesPval.csv")
DEG_monocle3_PDgene_FDR_up <- data.frame(up_rare_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_rare_genes"], method="fdr"), 
	up_exp_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_exp_genes"], method="fdr"), 
	up_meth_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_meth_genes"], method="fdr"), 
	up_rareV_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_rareV_genes"], method="fdr"), 
	up_gwas_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_gwas_genes"], method="fdr"), 
	up_cnv_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_cnv_genes"], method="fdr"), 
	up_gene4PD_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_gene4PD_genes"], method="fdr"), 
	up_mala_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_mala_genes"], method="fdr"), 
	up_mala_high_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_mala_high_genes"], method="fdr"),
	up_malaALS_genes=p.adjust(DEG_monocle3_PDgene_pval_combined[,"up_malaALS_genes"], method="fdr"))
write.csv(DEG_monocle3_PDgene_FDR_up, file="./human_PD/PDgenes/deg.monocle3.PDgenesFDR.up.csv")


DEG_NOISeq_PDgene_pval <- list()
DEG_NOISeq_PDgene_raregenes <- list()
DEG_NOISeq_PDgene_expgenes <- list()
DEG_NOISeq_PDgene_methgenes <- list()
DEG_NOISeq_PDgene_rareVgenes <- list()
DEG_NOISeq_PDgene_gwasgenes <- list()
DEG_NOISeq_PDgene_cnvgenes <- list()
DEG_NOISeq_PDgene_gene4PDgenes <- list()
DEG_NOISeq_PDgene_malagenes <- list()
DEG_NOISeq_PDgene_malahighgenes <- list()
DEG_NOISeq_PDgene_malaALSgenes <- list()

for (subclass_label_id in RNA_subclass_forDEG) {
	print(paste0(subclass_label_id, " start!"))
	# read deg
	subclass_name <- gsub(" ", "_", subclass_label_id)
	deg_file <- paste0("./rds/RNA/diff/after_integra/NOISeq/subclass/", subclass_name, "_noiseqbio_deg.csv")
	celltype_age_terms_deg <- read.csv(deg_file, row.names=1)

	deg_list <- list()
	deg_list[["all"]] <- celltype_age_terms_deg %>% tibble::rownames_to_column("name") %>% pull(name) %>% toupper
	deg_list[["up"]] <- celltype_age_terms_deg %>% tibble::rownames_to_column("name") %>% dplyr::filter(log2FC<0) %>% pull(name) %>% toupper
	deg_list[["down"]] <- celltype_age_terms_deg %>% tibble::rownames_to_column("name") %>% dplyr::filter(log2FC>0) %>% pull(name) %>% toupper
	express_file <- paste0("./rds/RNA/diff/after_integra/pseudo_bulk/subclass/", subclass_name, "_pseudobulk_autochr.csv")
	express_genes <- expressed_gene(express_file, express=0, cell=1)
	pval_list <- c()
	for (type in c("all", "up", "down")) {
		for (genes in list(rare_genes, exp_genes, meth_genes, rareV_genes, gwas_genes, cnv_genes, gene4PD_genes, mala_genes, mala_high_genes, malaALS_genes)) {
			ovlp_len <- length(intersect(deg_list[[type]], genes))
			if (ovlp_len>0) {
				pval <- 1-phyper(length(intersect(deg_list[[type]], genes)), length(deg_list[[type]]), length(express_genes)-length(deg_list[[type]]), length(genes))
			} else {
				pval <- NA
			}			
			pval_list <- c(pval_list, pval)
		}
	}
	DEG_NOISeq_PDgene_pval[[subclass_label_id]] <- pval_list

	ovlp_raregenes <- intersect(deg_list[["up"]], rare_genes)
	ovlp_expgenes <- intersect(deg_list[["up"]], exp_genes)
	ovlp_methgenes <- intersect(deg_list[["up"]], meth_genes)
	ovlp_rareVgenes <- intersect(deg_list[["up"]], rareV_genes)
	ovlp_gwasgenes <- intersect(deg_list[["up"]], gwas_genes)
	ovlp_cnvgenes <- intersect(deg_list[["up"]], cnv_genes)
	ovlp_gene4PDgenes <- intersect(deg_list[["up"]], gene4PD_genes)
	ovlp_malagenes <- intersect(deg_list[["up"]], mala_genes)
	ovlp_malahighgenes <- intersect(deg_list[["up"]], mala_high_genes)
	ovlp_malaALSgenes <- intersect(deg_list[["up"]], malaALS_genes)
	DEG_NOISeq_PDgene_raregenes[[subclass_label_id]] <- ovlp_raregenes
	DEG_NOISeq_PDgene_expgenes[[subclass_label_id]] <- ovlp_expgenes
	DEG_NOISeq_PDgene_methgenes[[subclass_label_id]] <- ovlp_methgenes
	DEG_NOISeq_PDgene_rareVgenes[[subclass_label_id]] <- ovlp_rareVgenes
	DEG_NOISeq_PDgene_gwasgenes[[subclass_label_id]] <- ovlp_gwasgenes
	DEG_NOISeq_PDgene_cnvgenes[[subclass_label_id]] <- ovlp_cnvgenes
	DEG_NOISeq_PDgene_gene4PDgenes[[subclass_label_id]] <- ovlp_gene4PDgenes
	DEG_NOISeq_PDgene_malagenes[[subclass_label_id]] <- ovlp_malagenes
	DEG_NOISeq_PDgene_malahighgenes[[subclass_label_id]] <- ovlp_malahighgenes
	DEG_NOISeq_PDgene_malaALSgenes[[subclass_label_id]] <- ovlp_malaALSgenes
}
DEG_NOISeq_PDgene_pval_combined <- do.call(rbind, DEG_NOISeq_PDgene_pval)
colnames(DEG_NOISeq_PDgene_pval_combined) <- c("all_rare_genes", "all_exp_genes", "all_meth_genes", "all_rareV_genes", "all_gwas_genes", "all_cnv_genes", "all_gene4PD_genes", "all_mala_genes", "all_mala_high_genes", "all_malaALS_genes", 
	"up_rare_genes", "up_exp_genes", "up_meth_genes", "up_rareV_genes", "up_gwas_genes", "up_cnv_genes", "up_gene4PD_genes", "up_mala_genes", "up_mala_high_genes", "up_malaALS_genes",  
	"down_rare_genes", "down_exp_genes", "down_meth_genes", "down_rareV_genes", "down_gwas_genes", "down_cnv_genes", "down_gene4PD_genes", "down_mala_genes", "down_mala_high_genes", "down_malaALS_genes")
write.csv(DEG_NOISeq_PDgene_pval_combined, file="./human_PD/PDgenes/deg.NOISeq.PDgenesPval.csv")
DEG_NOISeq_PDgene_FDR_up <- data.frame(up_rare_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_rare_genes"], method="fdr"), 
	up_exp_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_exp_genes"], method="fdr"), 
	up_meth_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_meth_genes"], method="fdr"), 
	up_rareV_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_rareV_genes"], method="fdr"), 
	up_gwas_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_gwas_genes"], method="fdr"), 
	up_cnv_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_cnv_genes"], method="fdr"), 
	up_gene4PD_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_gene4PD_genes"], method="fdr"), 
	up_mala_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_mala_genes"], method="fdr"), 
	up_mala_high_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_mala_high_genes"], method="fdr"),
	up_malaALS_genes=p.adjust(DEG_NOISeq_PDgene_pval_combined[,"up_malaALS_genes"], method="fdr"))
write.csv(DEG_NOISeq_PDgene_FDR_up, file="./human_PD/PDgenes/deg.NOISeq.PDgenesFDR.up.csv")

DEG_PDgene_FDR_up <- cbind(DEG_monocle3_PDgene_FDR_up[,c("up_gene4PD_genes", "up_mala_genes")], DEG_NOISeq_PDgene_FDR_up[,c("up_gene4PD_genes", "up_mala_genes")])


DEG_monocle3_PDgene_FDR_up %>% select("up_gene4PD_genes", "up_mala_genes") %>% tibble::rownames_to_column("subclass") -> DEG_monocle3_PDgene_FDR_select
DEG_NOISeq_PDgene_FDR_up %>% select("up_gene4PD_genes", "up_mala_genes") %>% tibble::rownames_to_column("subclass") -> DEG_NOISeq_PDgene_FDR_select

DEG_monocle3_PDgene_FDR_select %>% full_join(DEG_NOISeq_PDgene_FDR_select, by=join_by(subclass==subclass)) %>% arrange(up_gene4PD_genes.x, up_mala_genes.x, up_gene4PD_genes.y, up_mala_genes.y)

DEG_197 <- intersect(intersect(intersect(DEG_monocle3_PDgene_gene4PDgenes[["197 SNr Six3 Gaba"]], DEG_monocle3_PDgene_malagenes[["197 SNr Six3 Gaba"]]), DEG_NOISeq_PDgene_gene4PDgenes[["197 SNr Six3 Gaba"]]), DEG_NOISeq_PDgene_malagenes[["197 SNr Six3 Gaba"]])
DEG_327 <- intersect(intersect(intersect(DEG_monocle3_PDgene_gene4PDgenes[["327 Oligo NN"]], DEG_monocle3_PDgene_malagenes[["327 Oligo NN"]]), DEG_NOISeq_PDgene_gene4PDgenes[["327 Oligo NN"]]), DEG_NOISeq_PDgene_malagenes[["327 Oligo NN"]])
DEG_135 <- intersect(intersect(intersect(DEG_monocle3_PDgene_gene4PDgenes[["135 STN-PSTN Pitx2 Glut"]], DEG_monocle3_PDgene_malagenes[["135 STN-PSTN Pitx2 Glut"]]), DEG_NOISeq_PDgene_gene4PDgenes[["135 STN-PSTN Pitx2 Glut"]]), DEG_NOISeq_PDgene_malagenes[["135 STN-PSTN Pitx2 Glut"]])
DEG_215 <- intersect(intersect(intersect(DEG_monocle3_PDgene_gene4PDgenes[["215 SNc-VTA-RAmb Foxa1 Dopa"]], DEG_monocle3_PDgene_malagenes[["215 SNc-VTA-RAmb Foxa1 Dopa"]]), DEG_NOISeq_PDgene_gene4PDgenes[["215 SNc-VTA-RAmb Foxa1 Dopa"]]), DEG_NOISeq_PDgene_malagenes[["215 SNc-VTA-RAmb Foxa1 Dopa"]])
DEG_166 <- intersect(intersect(intersect(DEG_monocle3_PDgene_gene4PDgenes[["166 MRN Pou3f1 C1ql4 Glut"]], DEG_monocle3_PDgene_malagenes[["166 MRN Pou3f1 C1ql4 Glut"]]), DEG_NOISeq_PDgene_gene4PDgenes[["166 MRN Pou3f1 C1ql4 Glut"]]), DEG_NOISeq_PDgene_malagenes[["166 MRN Pou3f1 C1ql4 Glut"]])
DEG_139 <- intersect(intersect(intersect(DEG_monocle3_PDgene_gene4PDgenes[["139 PH-LHA Foxb1 Glut"]], DEG_monocle3_PDgene_malagenes[["139 PH-LHA Foxb1 Glut"]]), DEG_NOISeq_PDgene_gene4PDgenes[["139 PH-LHA Foxb1 Glut"]]), DEG_NOISeq_PDgene_malagenes[["139 PH-LHA Foxb1 Glut"]])
DEG_195 <- intersect(intersect(intersect(DEG_monocle3_PDgene_gene4PDgenes[["195 SNr-VTA Pax5 Npas1 Gaba"]], DEG_monocle3_PDgene_malagenes[["195 SNr-VTA Pax5 Npas1 Gaba"]]), DEG_NOISeq_PDgene_gene4PDgenes[["195 SNr-VTA Pax5 Npas1 Gaba"]]), DEG_NOISeq_PDgene_malagenes[["195 SNr-VTA Pax5 Npas1 Gaba"]])
DEG_334 <- intersect(intersect(intersect(DEG_monocle3_PDgene_gene4PDgenes[["334 Microglia NN"]], DEG_monocle3_PDgene_malagenes[["334 Microglia NN"]]), DEG_NOISeq_PDgene_gene4PDgenes[["334 Microglia NN"]]), DEG_NOISeq_PDgene_malagenes[["334 Microglia NN"]])
intersect(intersect(intersect(intersect(DEG_197,DEG_327),DEG_135),DEG_166),DEG_334)

# plot KEGG genes
expressed_gene <- function(file, express=0, cell=1) {
	express_mat <- read.csv(file, row.names=1)
	expressGene <- rownames(express_mat)[rowSums(express_mat>express)>=cell]
	return(expressGene)	
}

subclass_label_id <- "215 SNc-VTA-RAmb Foxa1 Dopa"
print(paste0(subclass_label_id, " start!"))
# read deg
subclass_name <- gsub(" ", "_", subclass_label_id)
deg_file <- paste0("./rds/RNA/diff/after_integra/NOISeq/subclass/", subclass_name, "_noiseqbio_deg.csv")
celltype_age_terms_deg <- read.csv(deg_file, row.names=1)

deg_up <- celltype_age_terms_deg %>% tibble::rownames_to_column("name") %>% dplyr::filter(log2FC<0) %>% pull(name)
express_file <- paste0("./rds/RNA/diff/after_integra/pseudo_bulk/subclass/", subclass_name, "_pseudobulk_autochr.csv")
express_genes <- expressed_gene(express_file, express=0, cell=1)
deg_GO_KEGG <- enrich_deg(deg_list, express_genes, keyType="ENTREZID", fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", ont="ALL", organism="mmu") 

toType_genes <- tryCatch({bitr(deg_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")}, error = function(e){NULL})
group_toTypes <- bitr(deg_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
deg_up <- group_toTypes[,"ENTREZID"]

background_enzids <- bitr(express_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
background_genes <- background_enzids[,"ENTREZID"]
groups_KEGGs <- enrichKEGG(deg_up, organism="mmu", universe=background_genes)
options(browser = "/home/kaw033/softwares/miniconda3/envs/scRNA/bin/firefox")
browseKEGG(groups_KEGGs, 'mmu05012')

library("pathview")
mmu05012 <- pathview(gene.data  = deg_up, pathway.id = "mmu05012", species = "mmu")
# Produce the native KEGG plot (PNG)
dme_mmu05012 <- pathview(gene.data=deg_up, pathway.id="mmu05012", species = "mmu")
# Produce a different plot (PDF) (not displayed here)
dmeF_mmu05012 <- pathview(gene.data=deg_up, pathway.id="mmu05012", species = "mmu", kegg.native = F)

detach(package:ClusterGVis)
detach(package:monocle)
library(monocle3)
cdsdata_celltype_deg_select <- cdsdata_celltype[rowData(cdsdata_celltype)$name %in% c("Fnbp1", "Apod", "Lrrc4c", "Cldn11"),]
pt_vln_deg_select <- plot_genes_violin(cdsdata_celltype_deg_select, group_cells_by="age", ncol=4) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(pt_vln_deg_select, file="./figures/fig2/monocle3_example_deg.pdf", width=20, height=4)





