suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
library("limma")
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
library("dendsort")
import::from(.from = "colors.R", .directory ="/package/R", SnapATACPalette)

# read data
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
Idents(RNA) <- RNA@meta.data$subclass_label_id
RNA$age <- factor(case_match(RNA$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18), levels=c(2, 6, 12, 18))
RNA$rep <- factor(case_match(RNA$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2"), levels=c("rep1", "rep2"))
RNA$age_rep <- factor(case_match(RNA$sampleID, "2m_rep1" ~ 2.1, "6m_rep1" ~ 6.1, "12m_rep1" ~ 12.1, "18m_rep1" ~ 18.1, "2m_rep2" ~ 2.2, "6m_rep2" ~ 6.2, "12m_rep2" ~ 12.2, "18m_rep2" ~ 18.2), levels=c(2.1,2.2,6.1,6.2,12.1,12.2,18.1,18.2))

KO_parse_filter <- readRDS("KEGG/KO_parse_filter.rds")

# filter sex chromosome
Autochr.genetable <- read.csv(file="rds/RNA/preprocess/Autochr.genetable.csv")
RNA <- subset(RNA, features=Autochr.genetable$gene_name)

RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)
RNA_subclass_forDEG <- intersect(rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell)>=100], rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0])

# pathway score
RNA <- AddModuleScore_UCell(RNA, features = KO_parse_filter)

# limma DEP
RNA_pathway_UCell <- RNA@meta.data[,c(paste0(names(KO_parse_filter), "_UCell"), "sampleID", "subclass_label_id", "age", "rep", "age_rep")]


DEP_limma <- function(RNA_pathway, subclass_label_id){
  RNA_pathway_subclass_label_id <- RNA_pathway[which(RNA_pathway$subclass_label_id==subclass_label_id),]
  design <- model.matrix(~as.numeric(as.character(RNA_pathway_subclass_label_id$age))+RNA_pathway_subclass_label_id$rep)
  ourData <- t(RNA_pathway_subclass_label_id[,1:length(KO_parse_filter)])
  fit <- lmFit(ourData, design)
  K <- rbind(c(0,1,0))
  cfit <- contrasts.fit(fit, t(K)) 
  efit <- eBayes(cfit, trend=TRUE)
  DEP_tbl <- topTable(efit, n=length(KO_parse_filter))
  return(DEP_tbl)
}

for (subclass_label_id in RNA_subclass_forDEG) {
	print(subclass_label_id)
	name <- gsub(" ", "_", subclass_label_id)
	DEP_tbl_UCell <- DEP_limma(RNA_pathway_UCell, subclass_label_id)
	write.csv(DEP_tbl_UCell, file=paste0("rds/RNA/diff/after_integra/pathway/UCell/Ucell_limma_", name, "_DEPathway.txt"), quote=FALSE)
	print(paste0(subclass_label_id, ": finished!"))
}

# heatmap
DEP_statistic <- function(file, var="adj.P.Val", rename=var){
	DEP_tbl <- read.csv(file)
	DEP_tbl_mod <- DEP_tbl %>% mutate(PathwayID=substr(X,1,8)) %>% arrange(PathwayID) %>% select_("PathwayID", var) 
	colnames(DEP_tbl_mod)[2] <- rename
	return(DEP_tbl_mod)
}

KO_filter <- readRDS("KEGG/KO_filter.rds")
KO_dict <- KO_filter %>% select_("L1_ID", "L1", "L2_ID", "L2", "L3_ID", "L3", "PathwayID") %>% distinct() %>% as.data.frame()
subclass_label_id.dict <- RNA@meta.data[,c("subclass_label_id", "class_label_id")] %>% distinct %>% filter(subclass_label_id %in% RNA_subclass_forDEG) %>% arrange(subclass_label_id)

DEP_dict_tbl <- KO_dict[match(sort(names(KO_parse_filter)), KO_dict$PathwayID), ]
rownames(DEP_dict_tbl) <- DEP_dict_tbl$PathwayID


cell_dict_tbl <- subclass_label_id.dict[match(colnames(DEP_tbl_modulescore_padj), subclass_label_id.dict$subclass_label_id), ]


# pathway score using UCell
DEP_tbl_UCell_padj <- data.frame(PathwayID=sort(names(KO_parse_filter)))
for (subclass_label_id in RNA_subclass_forDEG) {
	print(subclass_label_id)
	name <- gsub(" ", "_", subclass_label_id)
	file <- paste0("rds/RNA/diff/after_integra/pathway/UCell/Ucell_limma_", name, "_DEPathway.txt")
	DEP_tbl_UCell <- DEP_statistic(file, var="adj.P.Val", rename=subclass_label_id)
	DEP_tbl_UCell_padj <- plyr::join(DEP_tbl_UCell_padj, DEP_tbl_UCell, by="PathwayID")
}
DEP_tbl_UCell_padj <- DEP_tbl_UCell_padj %>% tibble::column_to_rownames("PathwayID") 
DEP_tbl_UCell_logpadj <- -log10(DEP_tbl_UCell_padj)


DEP_tbl_UCell_logFC <- data.frame(PathwayID=sort(names(KO_parse_filter)))
for (subclass_label_id in RNA_subclass_forDEG) {
	print(subclass_label_id)
	name <- gsub(" ", "_", subclass_label_id)
	file <- paste0("rds/RNA/diff/after_integra/pathway/UCell/Ucell_limma_", name, "_DEPathway.txt")
	DEP_tbl_UCell <- DEP_statistic(file, var="logFC", rename=subclass_label_id)
	DEP_tbl_UCell_logFC <- plyr::join(DEP_tbl_UCell_logFC, DEP_tbl_UCell, by="PathwayID")
}
DEP_tbl_UCell_logFC <- DEP_tbl_UCell_logFC %>% tibble::column_to_rownames("PathwayID") 

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
hclust_rows <- sort_hclust(hclust(dist(DEP_tbl_UCell_logFC), method = "ward.D2"))
hclust_cols <- sort_hclust(hclust(dist(t(DEP_tbl_UCell_logFC)),method = "ward.D2"))

DEP_labels <- tail(hclust_rows$order,20)
row_labels = structure(DEP_dict_tbl$L3, names=DEP_dict_tbl$PathwayID)
ha = rowAnnotation(L3=anno_mark(at=DEP_labels, labels=row_labels[DEP_labels]))

identical(colnames(DEP_tbl_UCell_logFC),colnames(DEP_tbl_UCell_padj))
identical(rownames(DEP_tbl_UCell_logFC),rownames(DEP_tbl_UCell_padj))

DEP_tbl_UCell_padj %>% mutate(across(everything(), ~case_when(.x < 0.001 ~ "***", .x < 0.01 ~ "**", .x < 0.05 ~ "*", .default = ""))) -> DEP_tbl_UCell_star


pt_DEP_UCell_logFCheat_cluster_label_dendrev_order <- Heatmap(DEP_tbl_UCell_logFC, name="Beta coefficient", 
#	col=col_fun, 
	column_order=c(1:ncol(DEP_tbl_UCell_logFC)), 
	#row_order=rows_order,
	cluster_rows=rev(as.dendrogram(hclust_rows)), 
	#cluster_columns=hclust_cols,
    show_row_names=TRUE, 
    row_names_side = "left",
    row_names_gp=gpar(fontsize=4),
    show_column_names=TRUE,
    column_names_gp=gpar(fontsize=8),
    top_annotation=ha_col,
    right_annotation = ha,
    left_annotation=ha_row,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(DEP_tbl_UCell_star[i, j], x, y, gp = gpar(fontsize = 5))}
)

pdf("./rds/RNA/diff/after_integra/pathway/UCell/UCell_DEP.logFC_padj.heatmap.cluster.label.pdf", width=12, height=18)
pt_DEP_UCell_logFCheat_cluster_label_dendrev_order
dev.off()

# statistic
PD_KEGG <- KO_parse_filter$mmu05012


