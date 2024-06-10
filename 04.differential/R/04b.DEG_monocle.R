#!/usr/bin/env Rscript

# read 
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
RNA$age <- case_match(RNA$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18)
RNA$rep <- case_match(RNA$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2")


# differential expression genes
expression_matrix <- RNA@assays$RNA@counts
cells <- RNA@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "integrated_snn_res.0.2", "sampleID", "class_label_id", "subclass_label_id", "age", "rep")]
cells$barcode <- str_split(rownames(cells), "_", simplify=TRUE)[,1]
write.csv(cells, file="./rds/RNA/diff/after_integra/monocle3/cell2anno.csv", quote=FALSE, row.names=FALSE)
genes <- data.frame(name=RNA@assays$RNA@counts@Dimnames[[1]], gene_short_name=RNA@assays$RNA@counts@Dimnames[[1]])
rownames(genes) <- RNA@assays$RNA@counts@Dimnames[[1]]
# Make the CDS object
RNA.cds <- new_cell_data_set(expression_matrix,
                         cell_metadata=cells,
                         gene_metadata=genes)
saveRDS(RNA.cds, file="./rds/RNA/diff/after_integra/monocle3/RNA.cds.rds")

Autochr.genetable <- read.csv(file="rds/RNA/preprocess/Autochr.genetable.csv")
RNA.cds.autochr <- RNA.cds[rowData(RNA.cds)$name %in% Autochr.genetable$gene_name,]

celltype_statistic <- function(cdsdata, column, celltype){
	cdsdata_celltype <- cdsdata[,colData(cdsdata)[,column] == celltype]
	celltype_acc <- Matrix::rowSums(exprs(cdsdata_celltype)>0)
	celltype_expressed_genes <- names(celltype_acc)[which(celltype_acc>10)]
	cdsdata_celltype <- cdsdata_celltype[rowData(cdsdata_celltype)$name %in% celltype_expressed_genes,]
	celltype_depth <- Matrix::colSums(exprs(cdsdata_celltype))
	cdsdata_celltype[["log10dep"]] <- log10(celltype_depth)

	gene_number <- nrow(cdsdata_celltype)
	cell_number <- ncol(cdsdata_celltype)
	cell_table <- table(cdsdata_celltype$age, cdsdata_celltype$rep)
	cell_samples <- as.numeric(cell_table)
	celltype_statistics <- c(celltype, gene_number, cell_number, cell_samples)
	names(celltype_statistics) <- c("celltype", "gene_number", "cell_number",
		paste(rep(rownames(cell_table), ncol(cell_table)), rep(colnames(cell_table), each=nrow(cell_table)), sep="_"))
	return(celltype_statistics)
}

subclass_gene_cell_stat <- c("subclass_label_id", "gene_number", "cell_number", paste(rep(c(2,6,12,18),2), rep(c("rep1", "rep2"), each=4), sep="_"))

for (i in levels(RNA.cds.autochr$subclass_label_id)) {
	subclass_gene_cell <- celltype_statistic(RNA.cds.autochr, column="subclass_label_id", celltype=i)
	subclass_gene_cell_stat <- rbind(subclass_gene_cell_stat, subclass_gene_cell)
}
subclass_gene_cell_stat <- subclass_gene_cell_stat[-1,]
rownames(subclass_gene_cell_stat) <- subclass_gene_cell_stat[,1]
subclass_gene_cell_stat <- as.data.frame(subclass_gene_cell_stat[,-1])

write.csv(subclass_gene_cell_stat, file="rds/RNA/diff/after_integra/monocle3/subclass_gene_cell_stat.csv")

RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)
RNA_subclass_forDEG <- intersect(rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell)>=100], rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0])

DEG_fun <- function(cdsdata, column, celltype, output){
	cdsdata_celltype <- cdsdata[,colData(cdsdata)[,column] == celltype]
	celltype_acc <- Matrix::rowSums(exprs(cdsdata_celltype)>0)
	celltype_expressed_genes <- names(celltype_acc)[which(celltype_acc>10)]
	cdsdata_celltype <- cdsdata_celltype[rowData(cdsdata_celltype)$name %in% celltype_expressed_genes,]
	celltype_depth <- Matrix::colSums(exprs(cdsdata_celltype))
	cdsdata_celltype[["log10dep"]] <- log10(celltype_depth)

	celltype_gene_fits <- fit_models(cdsdata_celltype, model_formula_str = "~age + log10dep + rep")
	celltype_fit_coefs <- coefficient_table(celltype_gene_fits)
	celltype_age_terms <- celltype_fit_coefs %>% filter(term == "age")
	celltype_age_terms_deg_allgene <- celltype_age_terms %>% select(name, num_cells_expressed, term, p_value, q_value, estimate)
	celltype_age_terms_deg <- celltype_age_terms %>% filter (q_value < 0.05) %>% select(name, num_cells_expressed, term, p_value, q_value, estimate) 
	cdsdata_celltype_deg <- cdsdata_celltype[rowData(cdsdata_celltype)$name %in% celltype_age_terms_deg$name,]
	celltype_rename <- gsub(" ", "_", celltype)
	saveRDS(celltype_gene_fits, file=paste0(output, celltype_rename, ".gene_fits.rds"))
	write.csv(celltype_age_terms_deg_allgene, file=paste0(output, celltype_rename, ".age_terms_allgenes.csv"))
	write.csv(celltype_age_terms_deg, file=paste0(output, celltype_rename, ".age_terms_deg.csv"))
	if (nrow(cdsdata_celltype_deg)==0) {
		message(celltype_rename, " has no DEG!")
	} else if (nrow(cdsdata_celltype_deg)>100) {
		cdsdata_celltype_deg <- cdsdata_celltype_deg[1:100,]
		pt_vln_deg <- plot_genes_violin(cdsdata_celltype_deg, group_cells_by="age", ncol=5) + theme(axis.text.x=element_text(angle=45, hjust=1))
		ggsave(pt_vln_deg, file=paste0(output, celltype_rename, ".deg_vln.pdf"), width=24, height=ceiling(nrow(cdsdata_celltype_deg)/5)+1)
	} else {
		pt_vln_deg <- plot_genes_violin(cdsdata_celltype_deg, group_cells_by="age", ncol=5) + theme(axis.text.x=element_text(angle=45, hjust=1))
		ggsave(pt_vln_deg, file=paste0(output, celltype_rename, ".deg_vln.pdf"), width=24, height=ceiling(nrow(cdsdata_celltype_deg)/5)+1)			
	}
}

for (i in RNA_subclass_forDEG) {
	print(i)
	DEG_fun(RNA.cds.autochr, column="subclass_label_id", celltype=i, output="rds/RNA/diff/after_integra/monocle3/subclass/")
}




