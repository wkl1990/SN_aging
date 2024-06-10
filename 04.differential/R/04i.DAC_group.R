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
library("ChIPseeker")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
options('homer_path' = "/home/kaw033/softwares/homer")
library("marge")


# read data 
input <- "rds/ATAC/after_integra/seurat/SN_integra.final.peak.srt.pmat.seurat4.rds"
ATAC <- readRDS(input)
Idents(ATAC) <- ATAC$subclass_label_id
ATAC$sampleID <- factor(ATAC$orig.ident, levels=c("2m_rep1", "2m_rep2", "6m_rep1", "6m_rep2", "12m_rep1", "12m_rep2", "18m_rep1", "18m_rep2"))
ATAC$age <- case_match(ATAC$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18)
ATAC$rep <- case_match(ATAC$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2")

ATAC_subclass_sample_cell <- table(ATAC$subclass_label_id, ATAC$sampleID)
ATAC_subclass_forDEG <- intersect(rownames(ATAC_subclass_sample_cell)[rowSums(ATAC_subclass_sample_cell)>=100], rownames(ATAC_subclass_sample_cell)[rowSums(ATAC_subclass_sample_cell<10)==0])

GO_database <- 'org.Mm.eg.db'
KEGG_database <- 'mmu'

# pesudobulk matrix
pseudo_feature <- function(data, ident=NULL, group_by="age", assay="ATAC", slot="count", deg_file, gene_column="name") {
	if (is.null(ident)) {
		data_subset <- data
	} else {
		data_subset <- subset(x=data, idents=ident)
	}
	bulk_data <- AggregateExpression(data_subset, group.by=group_by, assays=assay, slot=slot)
	bulk_data_ATAC <- as.matrix(bulk_data$ATAC)
	bulk_data_cpm <- edgeR::cpm(bulk_data_ATAC)
	deg_stat <- read.csv(deg_file, row.names=1)
	if (!(gene_column %in% colnames(deg_stat))) {
		deg_stat %>% tibble::rownames_to_column(gene_column) -> deg_stat
	}
	bulk_data_ATAC <- bulk_data_cpm[which(rownames(bulk_data_cpm) %in% deg_stat[,gene_column]),]
	return(bulk_data_ATAC)	
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


peak2gene <- function(peaks, tss=c(-1000,1000), flank=3000, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene) {
	peaks %>% as.data.frame() %>% dplyr::rename(peaks=".") %>% tidyr::separate(peaks, c("chr", "start", "end"), sep="[:-]") -> peaks_df
	peaks_GR <- makeGRangesFromDataFrame(peaks_df)
	peaks_gene <- tryCatch(seq2gene(peaks_GR, tssRegion=tss, flankDistance=flank, TxDb=txdb), error=function(e) NA)
	return(peaks_gene)
}

expressed_gene <- function(file, express=0, cell=1) {
	express_mat <- read.csv(file, row.names=1)
	expressGene <- rownames(express_mat)[rowSums(express_mat>express)>=cell]
	return(expressGene)	
}

peak_expressed_gene <- function(peak_genes, expressed_genes, orgdb=org.Mm.eg.db){
	expressed_genes_entrezid <- bitr(expressed_genes, "SYMBOL", "ENTREZID", OrgDb=orgdb)
	peak_genes_express <- peak_genes[which(peak_genes %in% expressed_genes_entrezid$ENTREZID)]
	return(peak_genes_express)
}

get_cluster_enzids <- function(cm, expressed_genes, orgdb=org.Mm.eg.db, flank=3000) {
	cm_cluster_enzids <- list()
	for (cluster in unique(cm$wide.res$cluster)) {
		cluster_name <- paste0("C", cluster)
		cm_cluster_peaks <- cm$wide.res$gene[which(cm$wide.res$cluster==cluster)]
		cm_cluster_peaks_genes <- peak2gene(cm_cluster_peaks, flank=flank)
		cm_cluster_peaks_express_genes <- peak_expressed_gene(cm_cluster_peaks_genes, expressed_genes, orgdb=orgdb)
		cm_cluster_enzids[[cluster_name]] <- cm_cluster_peaks_express_genes
	}
	return(cm_cluster_enzids)	
}

enrich_group <- function(group_genes, background_genes, keyType="ENTREZID", fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db", ont="ALL", organism="mmu") {
	groups <- 1:length(group_genes)
	group_names <- names(group_genes)
	if (keyType==fromType) {
		print("Attention: keyType is same as fromType, no conversion is needed!")
	} else if (keyType==toType) {
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

peak2df <- function(peaks) {
	peaks %>% as.data.frame() %>% dplyr::rename(peaks=".") %>% tidyr::separate(peaks, c("chr", "start", "end"), sep="[:-]") -> peaks_df
	return(peaks_df)
}

enrich_group_motif <- function(group_peaks, genome="mm10", size="given", core=30, homer_results_path) {
	groups <- 1:length(group_peaks)
	group_names <- names(group_peaks)
	group_dfs <-list()
	for (group_name in group_names) {
		group_dfs[[group_name]] <- peak2df(group_peaks[[group_name]])
	}
	groups_motifs <-list()
	for (group_name in group_names) {
		find_motifs_genome(group_dfs[[group_name]], path = file.path(homer_results_path, group_name), genome = genome, scan_size = size, optimize_count = 25, cores = core)
		groups_motifs[[group_name]]  <- read_known_results(path = file.path(homer_results_path, group_name), homer_dir = TRUE) %>% as.data.frame() %>% select(-motif_pwm)
		if (nrow(groups_motifs[[group_name]]) > 0) {
			groups_motifs[[group_name]]$group <- group_name
		}
	}
	group_motif <- do.call(rbind, groups_motifs)
	group_motif <- group_motif %>% mutate(pvalue=10**(-log_p_value)) %>% group_by(group) %>% mutate(rank=dense_rank(pvalue)) 
	return(group_motif)
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

group_plot_motif <- function(cm, group_genes, group_annot, deg_order, topgene=1, topannot=5) {
	group_annot %>% select(group, motif_name, pvalue, tgt_pct) %>% rename(Description="motif_name", GeneRatio="tgt_pct")  %>% group_by(group) %>% arrange(pvalue) %>% slice(1:topannot) %>% filter(pvalue<0.05) -> enrich2
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

for (subclass_label_id in ATAC_subclass_forDEG) {
	print(paste0(subclass_label_id, " start!"))
	# read deg
	subclass_name <- gsub(" ", "_", subclass_label_id)
	deg_file <- paste0("./rds/ATAC/after_integra/diff/monocle3/", subclass_name, ".age_terms_deg.csv")
	celltype_age_terms_deg <- read.csv(deg_file, row.names=1)
	if (nrow(celltype_age_terms_deg)<10) {
		next
	}
	celltype_age_terms_deg %>% arrange(q_value) %>% pull(name) -> deg_order
	celltype_age_terms_deg_up <- celltype_age_terms_deg %>% dplyr::filter(estimate>0)
	celltype_age_terms_deg_down <- celltype_age_terms_deg %>% dplyr::filter(estimate<0)
	# pesudobulk matrix (cpm and avg)
	ATAC_subset <- subset(x=ATAC, idents=subclass_label_id)
	ATAC_subset_bulk_ATAC_cpm <- pseudo_feature(ATAC, ident=subclass_label_id, deg_file=deg_file)
	# clustering
	# using mfuzz for clustering (mfuzz)
	cm_cpm <- clusterData(exp=ATAC_subset_bulk_ATAC_cpm, cluster.method="mfuzz", cluster.num=pre_cluster)
	cluster_genes_list <- get_cluster_genes(cm_cpm)
	group_cm_gene <- cm_cpm$wide.res %>% select(gene, cluster) %>% mutate(deg=case_when(gene %in% celltype_age_terms_deg_up$name ~ "up", gene %in% celltype_age_terms_deg_down$name ~ "down", .default="NA"))
	group_deg <- table(group_cm_gene$cluster,group_cm_gene$deg)
	write.csv(group_deg, file=paste0("./rds/ATAC/after_integra/diff/monocle3/group/group", pre_cluster, "/", subclass_name, ".deg2groups.csv"))
	express_file <- paste0("./rds/RNA/diff/after_integra/pseudo_bulk/subclass/", subclass_name, "_pseudobulk_autochr.csv")
	express_genes <- expressed_gene(express_file, express=0, cell=1)
	express_genes_enzid <- bitr(express_genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
	cluster_enzids_list <- get_cluster_enzids(cm_cpm, express_genes, orgdb=org.Mm.eg.db, flank=3000)
	group_GO_KEGG_list <- enrich_group(cluster_enzids_list, express_genes_enzid$ENTREZID, keyType="ENTREZID", fromType="ENTREZID", toType="None", OrgDb="org.Mm.eg.db", ont="ALL", organism="mmu") 
	write.csv(group_GO_KEGG_list, file=paste0("./rds/ATAC/after_integra/diff/monocle3/group/group", pre_cluster, "/", subclass_name, "_", cluster, "_group_GO_KEGG.csv")) 
	group_motif_list <- enrich_group_motif(cluster_genes_list, genome="mm10", size="given", core=30, homer_results_path=paste0("./rds/ATAC/after_integra/diff/monocle3/group/group", pre_cluster, "/homer_motif/", subclass_name, "_", cluster))
	write.csv(group_motif_list, file=paste0("./rds/ATAC/after_integra/diff/monocle3/group/group", pre_cluster, "/", subclass_name, "_", cluster, "_group_motif.csv")) 
	pt_visC_list <- NULL
	pt_visC_motif_list <- NULL
	if (nrow(group_GO_KEGG_list)>0) {
		pt_visC_list <- group_plot(cm_cpm, cluster_genes_list, group_GO_KEGG_list, deg_order, topgene=1, topannot=3)
	}
	if (nrow(group_motif_list)>0) {
		pt_visC_motif_list <- group_plot_motif(cm, cluster_genes_list, group_motif_list, deg_order, topgene=1, topannot=5)
	}
	pdf(paste0("./rds/ATAC/after_integra/diff/monocle3/group/group", pre_cluster, "/", subclass_name, ".test.groups.pdf"), width=16, height=9)
	if (!is.null(pt_visC_list)) {
		print(pt_visC_list)
	}
	if (!is.null(pt_visC_motif_list)) {
		print(pt_visC_motif_list)
	}
	dev.off()
	print(paste0(subclass_label_id, " finished!"))
#	rm()
}
