#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-r", "--resolution", default=NULL, help="resolution for clustering")
parser$add_argument("-a", "--assay", default="integrated", help="assay")
parser$add_argument("-m", "--method", default=4, help="clustering method, default Leiden")
parser$add_argument("-t", "--plot", action='store_true', help="whether to output table and plot")
parser$add_argument("-s", "--sample", default="sampleID", help="sample name in meta table")
parser$add_argument("-g", "--gene", default="Ptbp1,Ptbp2", help="gene name to plot")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
if (!is.null(args$resolution)){
  resolution = as.character(args$resolution)
}
assay = as.character(args$assay)
method = as.numeric(args$method)
uplot = ifelse(args$plot, TRUE, FALSE) 
sampleID = as.character(args$sample)
gene = as.character(args$gene)
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))

cluster_rna <- function(RNAdata, resolution=resolution, algorithm=method, assay="integrated") {
  DefaultAssay(RNAdata) <- assay
  RNAdata <- FindClusters(object=RNAdata, resolution=resolution, algorithm=algorithm)
  return(RNAdata)
}


print("Step1: read rds data!")
RNAdata <- readRDS(input)
if (exists(quote(resolution))) {
  resolution <- as.numeric(unlist(str_split(resolution, ",")))
} else {
  print("Resolution(s) not provided, use 1.0 default resolution!")
  resolution <- rep(1.0, length(RNAdata))
}
if (length(resolution)>1) {
  names(resolution) <- names(RNAdata)
}

print("Step2: run clustering!")
if (length(RNAdata)==1) {
  RNAdata_clusterlist <- cluster_rna(RNAdata, resolution=resolution, algorithm=method, assay=assay)
  RNAsceList <- RNAdata_clusterlist
} else {
  RNAsceList <- list()
  for (i in names(RNAdata)) {
    RNAdata_clusterlist <- cluster_rna(RNAdata[[i]], resolution=resolution[[i]], algorithm=method, assay=assay)
    RNAsceList[[i]] <- RNAdata_clusterlist
    rm(RNAdata_clusterlist)
  }
}

#Visualization
if (uplot) {
  print("Step3: output cluster table and umap plot!")
  if (length(RNAdata)==1) {
    RNA.clusters <- levels(Idents(RNAsceList))
    cluster.table <- table(Idents(RNAsceList))
    write.csv(cluster.table, file=paste0(outF, "_cluster.table.csv"))
    if (is.null(RNAsceList@meta.data[[sampleID]])) {
      RNAsceList@meta.data[[sampleID]] <- factor(RNAsceList@meta.data$orig.ident, levels=c("2m_rep1", "2m_rep2", "6m_rep1", "6m_rep2", "12m_rep1", "12m_rep2", "18m_rep1", "18m_rep2"))
    }
    pt_sample <- DimPlot(RNAsceList, reduction="umap", group.by=sampleID)
    ggsave(paste0(outF, ".", sampleID, ".umap.pdf"), plot=pt_sample, width=8, height=8)
    pt_splitdis <- DimPlot(RNAsceList, reduction="umap", split.by=sampleID, label=TRUE, label.size=4, ncol=2) + NoLegend()
    ggsave(paste0(outF, ".", sampleID, ".split.umap.pdf"), plot=pt_splitdis, width=8, height=4 * ceiling(length(unique(RNAsceList@meta.data[[sampleID]])) / 2))
    p_umap0 <- DimPlot(object=RNAsceList, reduction="umap", label=TRUE) + NoLegend()
    ggsave(p_umap0, file=paste0(outF, ".cluster.umap.pdf"), width=8, height=8)
    DefaultAssay(RNAsceList) <- "RNA"
    gene <- as.character(unlist(str_split(gene, ",")))
    pt_gene <- FeaturePlot(object=RNAsceList, features=gene, pt.size=0.1, max.cutoff="q95", reduction="umap", label=FALSE, ncol=2)
    ggsave(pt_gene, file=paste0(outF, ".gene.pdf"), width=8, height=4 * ceiling(length(gene) / 2))
  } else {
    for (i in names(RNAdata)) {
      RNA.clusters <- levels(Idents(RNAsceList[[i]]))
      cluster.table <- table(Idents(RNAsceList[[i]]))
      write.csv(cluster.table, file=paste0(outF, ".", i, "_cluster.table.csv"))

      if (is.null(RNAsceList[[i]]@meta.data[[sampleID]])) {
        RNAsceList[[i]]@meta.data[[sampleID]] <- factor(RNAsceList[[i]]@meta.data$orig.ident, levels=c("2m_rep1", "2m_rep2", "6m_rep1", "6m_rep2", "12m_rep1", "12m_rep2", "18m_rep1", "18m_rep2"))
      }
      pt_sample <- DimPlot(RNAsceList[[i]], reduction="umap", group.by=sampleID)
      ggsave(paste0(outF, ".", i, ".", sampleID, ".umap.pdf"), plot=pt_sample, width=8, height=8)
      pt_splitdis <- DimPlot(RNAsceList[[i]], reduction="umap", split.by=sampleID, label=TRUE, label.size=4, ncol=2) + NoLegend()
      ggsave(paste0(outF, ".", i, ".", sampleID, ".split.umap.pdf"), plot=pt_splitdis, width=8, height=4 * ceiling(length(unique(RNAsceList[[i]]@meta.data[[sampleID]])) / 2))
      p_umap0 <- DimPlot(object=RNAsceList[[i]], reduction="umap", label=TRUE) + NoLegend()
      ggsave(p_umap0, file=paste0(outF, ".", i, ".cluster.umap.pdf"), width=8, height=8)
      DefaultAssay(RNAsceList[[i]]) <- "RNA"
      pt_gene <- FeaturePlot(object=RNAsceList[[i]], features=gene, pt.size=0.1, max.cutoff="q95", reduction="umap", label=FALSE, ncol=3)
      ggsave(pt_gene, file=paste0(outF, ".", i, ".gene.pdf"), width=8, height=8 * ceiling(length(gene) / 3))
    }
  }
  print("Step4: save data!")
  saveRDS(RNAsceList, file=paste0(outF, ".cluster.rds"))
} else {
  print("Step3: save data!")
  saveRDS(RNAsceList, file=paste0(outF, ".cluster.rds"))
}
print("Job is done!")




