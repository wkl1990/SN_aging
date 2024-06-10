#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-c", "--column", default="celltype", help="celltype column")
parser$add_argument("-t", "--celltype", default=NULL, help="celltype to subset")
parser$add_argument("-p", "--cores", default=6, help="number of cores")
#parser$add_argument("-t", "--plot", action='store_true', help="whether to output table and plot")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
if (!is.null(args$celltype)){
  celltype = as.character(args$celltype)
}
column = as.character(args$column)
cores = as.numeric(args$cores)
#uplot = ifelse(args$plot, TRUE, FALSE) 
output = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("monocle3"))


#input <- "./rds/ATAC/preprocess/diff/monocle/ATAC.cds.rds"
ATAC.cds <- readRDS(input)
ATAC.cds.autochr <- ATAC.cds[grep("chr[1-9]", rowData(ATAC.cds)$name),]

cdsdata <- ATAC.cds.autochr
rm(ATAC.cds, ATAC.cds.autochr)
gc()
cdsdata_celltype <- cdsdata[,colData(cdsdata)[,column] == celltype]
celltype_acc <- Matrix::rowSums(exprs(cdsdata_celltype)>0)
celltype_expressed_genes <- names(celltype_acc)[which(celltype_acc>10)]
cdsdata_celltype <- cdsdata_celltype[rowData(cdsdata_celltype)$name %in% celltype_expressed_genes,]
celltype_depth <- Matrix::colSums(exprs(cdsdata_celltype))
cdsdata_celltype[["log10dep"]] <- log10(celltype_depth)

celltype_gene_fits <- fit_models(cdsdata_celltype, model_formula_str = "~age + log10dep + rep", cores=cores)
celltype_fit_coefs <- coefficient_table(celltype_gene_fits)
celltype_age_terms <- celltype_fit_coefs %>% filter(term == "age")
celltype_age_terms_deg <- celltype_age_terms %>% filter (q_value < 0.05) %>% select(name, num_cells_expressed, term, p_value, q_value, estimate) 
cdsdata_celltype_deg <- cdsdata_celltype[rowData(cdsdata_celltype)$name %in% celltype_age_terms_deg$name,]
if (nrow(cdsdata_celltype_deg)>100) {
	cdsdata_celltype_deg <- cdsdata_celltype_deg[1:100,]
}
pt_vln_deg <- plot_genes_violin(cdsdata_celltype_deg, group_cells_by="age", ncol=5) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(pt_vln_deg, file=paste0(output, gsub("/", "_", celltype), ".deg_vln.pdf"), width=24, height=ceiling(nrow(cdsdata_celltype_deg)/5)+1)
write.csv(celltype_age_terms_deg, file=paste0(output, gsub("/", "_", celltype), ".age_terms_deg.csv"))
saveRDS(celltype_gene_fits, file=paste0(output, gsub("/", "_", celltype), ".gene_fits.rds"))


