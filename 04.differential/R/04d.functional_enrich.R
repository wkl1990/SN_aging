suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("monocle3"))

# read data
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
RNA$age <- case_match(RNA$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18)
RNA$rep <- case_match(RNA$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2")

RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)
RNA_subclass_forDEG <- intersect(rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell)>=100], rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0])


# GO and KEGG for monocle3 deg
enrich_deg <- function(deg_genes, background_genes, keyType="ENTREZID", fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", ont="ALL", organism="mmu") {
	groups <- 1:length(deg_genes)
	group_names <- names(deg_genes)
	if (keyType==toType) {
		group_toTypes <-list()
		for (group_name in group_names) {
			toType_genes <- tryCatch({bitr(deg_genes[[group_name]], fromType=fromType, toType=toType, OrgDb=OrgDb)}, error = function(e){NULL})
			if (is.null(toType_genes)) {next}
			group_toTypes[[group_name]] <- bitr(deg_genes[[group_name]], fromType=fromType, toType=toType, OrgDb=OrgDb)
			deg_genes[[group_name]] <- group_toTypes[[group_name]][,toType]
		}
		background_enzids <- bitr(background_genes, fromType=fromType, toType=toType, OrgDb=OrgDb)
		background_genes <- background_enzids[,toType]
	} else if (keyType!=fromType) {
		stop("Error: keyType should be same as fromType or toType!")
	}
	groups_GOs <-list()
	for (group_name in group_names) {
		groups_GOs[[group_name]] <- enrichGO(deg_genes[[group_name]], OrgDb=OrgDb, keyType=keyType, ont=ont, universe=background_genes)
		groups_GOs[[group_name]] <- as.data.frame(groups_GOs[[group_name]])
		if (nrow(groups_GOs[[group_name]]) > 0) {
			groups_GOs[[group_name]]$group <- group_name
			groups_GOs[[group_name]]$DB <- "GO"
		}
	}
	group_GO <- do.call(rbind, groups_GOs)
	if (nrow(group_GO)>0) {group_GO <- group_GO %>% rename(subcategory = ONTOLOGY)}
	# KEGG can only use ENTREZID
	groups_KEGGs <-list()
	for (group_name in group_names) {
		groups_KEGGs[[group_name]] <- enrichKEGG(deg_genes[[group_name]], organism=organism, universe=background_genes)
		groups_KEGGs[[group_name]] <- as.data.frame(groups_KEGGs[[group_name]])
		if (nrow(groups_KEGGs[[group_name]]) > 0) {
			groups_KEGGs[[group_name]]$group <- group_name
			groups_KEGGs[[group_name]]$DB <- "KEGG"
		}
	}
	group_KEGG <- do.call(rbind, groups_KEGGs)
	group_KEGG$category <- NULL
	group_GO_KEGG <- rbind(group_GO, group_KEGG)
	return(group_GO_KEGG)
}

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
	deg_list[["all"]] <- celltype_age_terms_deg %>% pull(name)
	deg_list[["up"]] <- celltype_age_terms_deg %>% dplyr::filter(estimate>0) %>% pull(name)
	deg_list[["down"]] <- celltype_age_terms_deg %>% dplyr::filter(estimate<0) %>% pull(name)
	express_file <- paste0("./rds/RNA/diff/after_integra/pseudo_bulk/subclass/", subclass_name, "_pseudobulk_autochr.csv")
	express_genes <- expressed_gene(express_file, express=0, cell=1)
	deg_GO_KEGG <- enrich_deg(deg_list, express_genes, keyType="ENTREZID", fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", ont="ALL", organism="mmu") 
	write.csv(deg_GO_KEGG, file=paste0("./rds/RNA/diff/after_integra/monocle3/subclass/", subclass_name, ".deg.GO_KEGG.csv"))
}

# GO and KEGG for NOISeq deg
for (subclass_label_id in RNA_subclass_forDEG) {
	print(paste0(subclass_label_id, " start!"))
	# read deg
	subclass_name <- gsub(" ", "_", subclass_label_id)
	deg_file <- paste0("./rds/RNA/diff/after_integra/NOISeq/subclass/", subclass_name, "_noiseqbio_deg.csv")
	celltype_age_terms_deg <- read.csv(deg_file, row.names=1)
#	if (nrow(celltype_age_terms_deg)<10) {
#		next
#	}
	deg_list <- list()
	deg_list[["all"]] <- celltype_age_terms_deg %>% tibble::rownames_to_column("name") %>% pull(name)
	deg_list[["up"]] <- celltype_age_terms_deg %>% tibble::rownames_to_column("name") %>% dplyr::filter(log2FC<0) %>% pull(name)
	deg_list[["down"]] <- celltype_age_terms_deg %>% tibble::rownames_to_column("name") %>% dplyr::filter(log2FC>0) %>% pull(name)
	express_file <- paste0("./rds/RNA/diff/after_integra/pseudo_bulk/subclass/", subclass_name, "_pseudobulk_autochr.csv")
	express_genes <- expressed_gene(express_file, express=0, cell=1)
	deg_GO_KEGG <- enrich_deg(deg_list, express_genes, keyType="ENTREZID", fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", ont="ALL", organism="mmu") 
	write.csv(deg_GO_KEGG, file=paste0("./rds/RNA/diff/after_integra/NOISeq/subclass/", subclass_name, ".deg.GO_KEGG.csv"))
}



