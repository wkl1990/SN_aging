suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
library("ClusterGVis")
library("org.Mm.eg.db")
library("clusterProfiler")


# read data 
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
RNA$age <- case_match(RNA$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18)
RNA$rep <- case_match(RNA$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2")

RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)
RNA_subclass_forDEG <- intersect(rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell)>=100], rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0])


# check cluster
Idents(RNA) <- RNA$subclass_label_id
for (subclass_label_id in RNA_subclass_forDEG) {
	print(subclass_label_id)
	subclass_file <- gsub(" ", "_", subclass_label_id)
	subclass_age_terms_deg <- read.csv(paste0("./rds/RNA/diff/after_integra/monocle3/subclass/", subclass_file,".age_terms_deg.csv"), row.names=1)
	if (nrow(subclass_age_terms_deg)==0) {
		print(paste0(subclass_file, " has no deg!"))
		next
	}
	RNA_subset <- subset(x=RNA, idents=subclass_label_id)
	# pesudobulk matrix
	#RNA_subset_bulk <- AggregateExpression(RNA_subset, features=subclass_age_terms_deg$name, group.by="age")
	RNA_subset_bulk <- AverageExpression(RNA_subset, features=subclass_age_terms_deg$name, assays="RNA", group.by="age")
	RNA_subset_bulk_RNA <- as.matrix(RNA_subset_bulk$RNA)

	# check optimal cluster numbers
	tryCatch({
		pt_cluster <- getClusters(exp=RNA_subset_bulk_RNA)
		ggsave(pt_cluster, file=paste0("./rds/RNA/diff/after_integra/monocle3/group/subclass/", subclass_file, ".getCluster.pdf"), width=8, height=8)
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


# pesudobulk matrix
pseudo_feature <- function(data, ident=NULL, group_by="age", assay="RNA", slot="count", deg_file, gene_column="name") {
	if (is.null(ident)) {
		data_subset <- data
	} else {
		data_subset <- subset(x=data, idents=ident)
	}
	bulk_data <- AggregateExpression(data_subset, group.by=group_by, assays=assay, slot=slot)
	bulk_data_RNA <- as.matrix(bulk_data$RNA)
	bulk_data_cpm <- edgeR::cpm(bulk_data_RNA)
	deg_stat <- read.csv(deg_file, row.names=1)
	if (!(gene_column %in% colnames(deg_stat))) {
		deg_stat %>% tibble::rownames_to_column(gene_column) -> deg_stat
	}
	bulk_data_RNA <- bulk_data_cpm[which(rownames(bulk_data_cpm) %in% deg_stat[,gene_column]),]
	return(bulk_data_RNA)	
}

expressed_gene <- function(file, express=0, cell=1) {
	express_mat <- read.csv(file, row.names=1)
	expressGene <- rownames(express_mat)[rowSums(express_mat>express)>=cell]
	return(expressGene)	
}

get_cluster_genes <- function(cm) {
	cm_cluster_genes <- list()
	for (cluster in unique(cm$wide.res$cluster)) {
		cluster_name <- paste0("C", cluster)
		cm_cluster_genes[[cluster_name]] <- cm$wide.res$gene[which(cm$wide.res$cluster==cluster)]
	}
	return(cm_cluster_genes)	
}

Autochr.genetable <- read.csv(file="rds/RNA/preprocess/Autochr.genetable.csv")
Autochr.genetable_enzid <- bitr(Autochr.genetable$gene_name, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
write.csv(Autochr.genetable_enzid, file="rds/RNA/preprocess/Autochr.genetable.entrezid.csv")

enrich_group <- function(group_genes, background_genes, keyType="ENTREZID", fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", ont="ALL", organism="mmu") {
	groups <- 1:length(group_genes)
	group_names <- names(group_genes)
	if (keyType==toType) {
		group_toTypes <-list()
		for (group_name in group_names) {
			toType_genes <- tryCatch({bitr(group_genes[[group_name]], fromType=fromType, toType=toType, OrgDb=OrgDb)}, error = function(e){NULL})
			if (is.null(toType_genes)) {next}
			group_toTypes[[group_name]] <- bitr(group_genes[[group_name]], fromType=fromType, toType=toType, OrgDb=OrgDb)
			group_genes[[group_name]] <- group_toTypes[[group_name]][,toType]
		}
		background_enzids <- bitr(background_genes, fromType=fromType, toType=toType, OrgDb=OrgDb)
		background_genes <- background_enzids[,toType]
	} else if (keyType!=fromType) {
		stop("Error: keyType should be same as fromType or toType!")
	}
	groups_GOs <-list()
	for (group_name in group_names) {
		groups_GOs[[group_name]] <- enrichGO(group_genes[[group_name]], OrgDb=OrgDb, keyType=keyType, ont=ont, universe=background_genes)
		groups_GOs[[group_name]] <- as.data.frame(groups_GOs[[group_name]])
		if (nrow(groups_GOs[[group_name]]) > 0) {
			groups_GOs[[group_name]]$group <- group_name
			groups_GOs[[group_name]]$DB <- "GO"
		}
	}
	group_GO <- do.call(rbind, groups_GOs)
	group_GO <- group_GO %>% rename(subcategory = ONTOLOGY)
	# KEGG can only use ENTREZID
	groups_KEGGs <-list()
	for (group_name in group_names) {
		groups_KEGGs[[group_name]] <- enrichKEGG(group_genes[[group_name]], organism=organism, universe=background_genes)
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

group_plot <- function(cm, group_genes, group_annot, deg_order, topgene=1, topannot=3) {
	group_annot %>% arrange(group, DB, pvalue) %>% group_by(group, DB) %>% mutate(rank=dense_rank(pvalue), Description=gsub(" - Mus musculus \\(house mouse\\)", "", Description), Description=paste0(Description, " (", DB, ")")) %>% filter(rank %in% c(1:topannot)) %>% as.data.frame -> enrich2
	enrich2 %>% select(group, Description, pvalue, GeneRatio) -> enrich2
	groups <- 1:length(group_genes)
	group_names <- names(group_genes)
	markGenes <- c()
	for (group_name in group_names) {
		group_markGenes <- head(intersect(deg_order, group_genes[[group_name]]), topgene)
		markGenes <- c(markGenes, group_markGenes)
	}
	go_col=rep(ggsci::pal_d3()(length(unique(enrich2$group))),table(enrich2$group))
	pt_visC <- visCluster(object=cm, plot.type="both", ms.col=c("green","orange","red"), add.box=TRUE, boxcol=ggsci::pal_npg()(length(group_genes)), column_names_rot=45, sample.col=ggsci::pal_npg()(length(unique(cm$long.res$cell_type))), 
	  ctAnno.col=ggsci::pal_d3()(length(group_genes)), show_row_dend=F, markGenes=markGenes, markGenes.side="left", 
	  genes.gp=c('italic',fontsize=12,col="black"), annoTerm.data=enrich2, line.side="left", go.col=go_col)
	return(pt_visC)
}

#monocle3 deg
for (subclass_label_id in RNA_subclass_forDEG) {
	print(paste0(subclass_label_id, " start!"))
	# read deg
	subclass_name <- gsub(" ", "_", subclass_label_id)
	deg_file <- paste0("./rds/RNA/diff/after_integra/monocle3/subclass/", subclass_name, ".age_terms_deg.csv")
	celltype_age_terms_deg <- read.csv(deg_file, row.names=1)
	if (nrow(celltype_age_terms_deg)<10) {
		next
	}
	celltype_age_terms_deg %>% arrange(q_value) %>% pull(name) -> deg_order
	celltype_age_terms_deg_up <- celltype_age_terms_deg %>% dplyr::filter(estimate>0)
	celltype_age_terms_deg_down <- celltype_age_terms_deg %>% dplyr::filter(estimate<0)
	# pesudobulk matrix (cpm and avg)
	RNA_subset <- subset(x=RNA, idents=subclass_label_id)
	RNA_subset_bulk_RNA_cpm <- pseudo_feature(RNA, ident=subclass_label_id, deg_file=deg_file)
	# clustering
	# using mfuzz for clustering (mfuzz)
#	detach(package:monocle3)
	cm_cpm <- clusterData(exp=RNA_subset_bulk_RNA_cpm, cluster.method="mfuzz", cluster.num=pre_cluster)
	cluster_genes_list <- get_cluster_genes(cm_cpm)
	group_cm_gene <- cm_cpm$wide.res %>% select(gene, cluster) %>% mutate(deg=case_when(gene %in% celltype_age_terms_deg_up$name ~ "up", gene %in% celltype_age_terms_deg_down$name ~ "down", .default="NA"))
	group_deg <- table(group_cm_gene$cluster,group_cm_gene$deg)
	write.csv(group_deg, file=paste0("./rds/RNA/diff/after_integra/monocle3/group/subclass/group", pre_cluster, "/", subclass_name, ".deg2groups.csv"))
	express_file <- paste0("./rds/RNA/diff/after_integra/pseudo_bulk/subclass/", subclass_name, "_pseudobulk_autochr.csv")
	express_genes <- expressed_gene(express_file, express=0, cell=1)
	group_GO_KEGG_list <- enrich_group(cluster_genes_list, express_genes, keyType="ENTREZID", fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", ont="ALL", organism="mmu")
	write.csv(group_GO_KEGG_list, file=paste0("./rds/RNA/diff/after_integra/monocle3/group/subclass/group", pre_cluster, "/", subclass_name, "_group_GO_KEGG.csv")) 
	pt_visC_list <- group_plot(cm, cluster_genes_list, group_GO_KEGG_list, deg_order, topgene=1, topannot=3)
	pdf(paste0("./rds/RNA/diff/after_integra/monocle3/group/subclass/group", pre_cluster, "/", subclass_name, ".test.groups.pdf"), width=16, height=9)
	print(pt_visC_list)
	dev.off()
	print(paste0(subclass_label_id, " finished!"))
#	rm()
}


