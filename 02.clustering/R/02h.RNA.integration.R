#!/usr/bin/R

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-r", "--ref", required=TRUE, help="ref input RNA rds")
parser$add_argument("-q", "--query", required=TRUE, help="query input RNA rds")
parser$add_argument("-a1", "--ref_assay", required=TRUE, help="ref assay")
parser$add_argument("-a2", "--query_assay", required=TRUE, help="query assay")
parser$add_argument("-c1", "--ref_cluster", required=TRUE, help="ref cluster in meta table")
parser$add_argument("-c2", "--query_cluster", required=TRUE, help="query cluster in meta table")
parser$add_argument("-n", "--normalized", default="Strandard", help="normalization method, Strandard or SCTransform")
parser$add_argument("-i", "--integrated", default="cca", help="integrated method, cca or rpca")
parser$add_argument("-d", "--dynamic_downsample", action="store_true", help="dynamic downsampling number of cells for data")
parser$add_argument("-d1", "--ref_downsample", default=NULL, help="downsampling number of cells for ref data")
parser$add_argument("-d2", "--query_downsample", default=NULL, help="downsampling number of cells for query data")
parser$add_argument("-m", "--remove", default=NULL, help="filter celltype with small cell number in ref data")
parser$add_argument("-k", "--kanchor", default=5, help="number of anchors")
parser$add_argument("-f", "--nfeature", default=2000, help="number of selected feature")
parser$add_argument("-c", "--npc", default=25, help="number of pc used")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

#ref="/projects/ps-renlab2/kaw033/WK_aging/allen_data/2023Update/seurat/SN/allen_all_SN.seurat4.cluster.rds"; query="/projects/ps-renlab2/kaw033/WK_aging/rds/RNA/cluster/round1/RNA.combined.reidents.rds";
#ref_assay="RNA4"; query_assay="RNA"; ref_cluster="cellclass"; query_cluster="celltype"; normalized="Strandard"; integrated="cca"; ref_downsample=500; query_downsample=1000; remove_num=30; kanchor=30; nfeature=2000; npc=20;

ref = as.character(args$ref)
query = as.character(args$query)
ref_assay = as.character(args$ref_assay)
query_assay = as.character(args$query_assay)
ref_cluster = as.character(args$ref_cluster)
query_cluster = as.character(args$query_cluster)
normalized = as.character(args$normalized)
integrated = as.character(args$integrated)
dynamic_downsample = ifelse(args$dynamic_downsample, TRUE, FALSE)
if(!is.null(args$ref_downsample)) {ref_downsample = as.numeric(args$ref_downsample)}
if(!is.null(args$query_downsample)) {query_downsample = as.numeric(args$query_downsample)}
if(!is.null(args$remove)) {remove_num = as.numeric(args$remove)}
kanchor = as.numeric(args$kanchor)
nfeature = as.numeric(args$nfeature)
npc = as.numeric(args$npc)
outF = as.character(args$output)

#suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("circlize"))
library("dendsort")
#library("tictoc")

#packdir <- "./package/R"
#import::from(.from = "colors.R", .directory = packdir, ArchRPalettes)

# * load rds 
ref_data <- readRDS(ref)
query_data <- readRDS(query)
ref_data@active.assay <- ref_assay
query_data@active.assay <- query_assay

SetAssayNull <- function(data) {
  assays <- setdiff(names(data@assays), data@active.assay)
  if (length(assays) != 0) {
    for (assay in assays) {
      data[[assay]] <- NULL
    }
  }
  return(data)
} 
ref_data <- SetAssayNull(ref_data)
query_data <- SetAssayNull(query_data)

if(!(ref_cluster %in% names(ref_data@meta.data))) {
  ref_data@meta.data[,ref_cluster] <- Idents(ref_data)
}
if(!(query_cluster %in% names(query_data@meta.data))) {
  query_data@meta.data[,query_cluster] <- Idents(query_data)
}


removeSeurat <- function(seu, groupBy = "cluster_id", remove = 10) {
  groups <- unlist(seu[[groupBy]])
  groups_table <- table(groups)
  groups_rm <- names(groups_table)[groups_table<remove]
  cells <- colnames(seu)[!(groups %in% groups_rm)]
  #seu_rm <- subset(x=seu, subset=groupBy %in% groups_rm, invert=TRUE)
  seu_rm <- subset(x=seu, cells=cells)
  return(seu_rm)
}

if (exists("remove_num")) {
  ref_data <- removeSeurat(ref_data, groupBy = ref_cluster, remove = remove_num)
#  query_data <- removeSeurat(query_data, groupBy = query_cluster, remove = remove_num)
}


get.downsample.fun <- function(minNum = 500, maxNum = 1000) {
  f <- function(index, labels) {
    set.seed(2024)
    stat <- table(labels)
    ## put divided in the first since sometimes the full production may
    ## introduce overflow problem in R since
    # only 32-bit integer are supported
    idealNum <- (stat / sum(stat)) * minNum * length(stat)
    nms <- names(idealNum)
    # FIXME: logic is too complex
    lindex <- lapply(nms, function(i) {
      curIndex <- index[which(labels %in% i)]
      if ( (stat[i] <= max(minNum, idealNum[i])) ) {
        return(curIndex)
      }
      tmpIndex <- sample(curIndex, size = idealNum[i], replace = FALSE)
      tmp2Index <- tmpIndex
      if(length(tmpIndex) <= minNum) {
        tmp2Index <- sample(curIndex, size = minNum, replace = FALSE)
      }
      if (length(tmpIndex) > maxNum) {
        tmp2Index <- sample(tmpIndex, size = maxNum, replace = FALSE)
      }
      return(tmp2Index)
    })
    return(sort(unlist(lindex)))
  }
  return(f)
}

downsampleSeurat <- function(seu,
                             groupBy = "cluster_id",
                             minNum = 1000,
                             maxNum = 2000) {
  fn.dp <- get.downsample.fun(
    minNum = minNum, maxNum = maxNum)
  allcells <- colnames(seu)
  ## NOTE: seu[[groupBy]] will return data.frame with one column
  labels <- seu@meta.data[[groupBy]]
  if(is.null(labels)) {
    warning(groupBy, " is not int meta.data.")
    return(seu)
  }
  dp.cells <- fn.dp(index = allcells, labels = labels)
  subset(seu, cells = dp.cells)
}

if (dynamic_downsample) {
  ref_cellnum <- ncol(ref_data)
  query_cellnum <- ncol(query_data)
  ref_stat <- table(ref_data@meta.data[ref_cluster])
  query_stat <- table(query_data@meta.data[query_cluster])
  if (ref_cellnum<=query_cellnum) {
    ref_downsample <- round(median(ref_stat), 0)
    query_downsample <- round(ref_downsample*length(ref_stat)/length(query_stat),0)
  } else {
    query_downsample <- round(median(query_stat), 0)
    ref_downsample <- round(query_downsample*length(query_stat)/length(ref_stat),0)
  }
  message("Dynamic downsample: ", "ref_downsample: ", ref_downsample, ", query_downsample: ", query_downsample)
  message("Number of clusters: ", "ref clusters: ", length(unique(ref_data@meta.data[[ref_cluster]])), ", query clusters: ", length(unique(query_data@meta.data[[query_cluster]])))  
}

if (exists("ref_downsample")) {
  #message("require: downsample Seurat.")
  #message(str_glue("atac: {atac_dsMin} - {atac_dsMax}."))
  ref_data <- downsampleSeurat(seu = ref_data, groupBy = ref_cluster,
    minNum = ref_downsample, maxNum = ref_downsample)
}
if (exists("query_downsample")) {
  #message("require: downsample Seurat.")
  #message(str_glue("atac: {atac_dsMin} - {atac_dsMax}."))
  query_data <- downsampleSeurat(seu = query_data, groupBy = query_cluster,
    minNum = query_downsample, maxNum = query_downsample)
}


# * reference-based clustering
ref_data@meta.data$refCluster <- ref_data@meta.data[, ref_cluster]
query_data@meta.data$queryCluster <- query_data@meta.data[, query_cluster]

ref_data@meta.data$set <- "ref"
query_data@meta.data$set <- "query"

se.lst <- list(ref=ref_data, query=query_data)
rm(ref_data, query_data)
gc()

if (normalized=="Strandard") {
  # Strandard normalization
  # normalize and identify variable features for each dataset independently
  se.lst <- lapply(X=se.lst, FUN=function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method="vst", nfeatures=nfeature)
  })
  # select features that are repeatedly variable across datasets for integration
  se.features <- SelectIntegrationFeatures(object.list=se.lst, nfeatures=nfeature)
  se.anchors <- FindIntegrationAnchors(object.list=se.lst, anchor.features=se.features, k.anchor=kanchor, reduction=integrated)
  # this command creates an 'integrated' data assay
  se.combined <- IntegrateData(anchorset = se.anchors)
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(se.combined) <- "integrated"
  # Run the standard workflow for visualization and clustering
  se.combined <- ScaleData(se.combined, verbose = FALSE)
} else if(normalized=="SCTransform") {
  # SCTransform
  #se.lst <- lapply(X = se.lst, FUN = SCTransform)
  for (i in names(se.lst)) {
      se.lst[[i]] <- SCTransform(se.lst[[i]], verbose = TRUE, return.only.var.genes = FALSE, assay=se.lst[[i]]@active.assay)
  }
  se.features <- SelectIntegrationFeatures(object.list = se.lst, nfeatures = nfeature)
  se.lst <- PrepSCTIntegration(object.list = se.lst, anchor.features = se.features)
  se.anchors <- FindIntegrationAnchors(object.list = se.lst, normalization.method = "SCT", anchor.features = se.features, k.anchor=kanchor, reduction=integrated)
  se.combined <- IntegrateData(anchorset = se.anchors, normalization.method = "SCT")
} else {
    stop("Error: normalization method is not supported!")
}

se.combined <- RunPCA(se.combined, npcs = 100, verbose = FALSE)
pt_combined_elbow <- ElbowPlot(se.combined, ndims=100)
se.combined <- RunUMAP(se.combined, reduction = "pca", dims = 1:npc)
se.combined <- FindNeighbors(se.combined, reduction = "pca", dims = 1:npc)
se.combined <- FindClusters(se.combined, resolution = 0.1, algorithm = 1)
# Visualization
pt_coembed <- DimPlot(se.combined, reduction = "umap", group.by = "set")
pt_umap_refCluster <- DimPlot(se.combined, reduction = "umap", group.by = ref_cluster, label = TRUE, repel = FALSE) + NoLegend()
pt_umap_queryCluster <- DimPlot(se.combined, reduction = "umap", group.by = query_cluster, label = TRUE, repel = FALSE) + NoLegend()
se.combined$Cluster <- ifelse(is.na(se.combined$refCluster), se.combined$queryCluster, se.combined$refCluster)
pt_umap_combineCluster <- DimPlot(se.combined, reduction = "umap", group.by = "Cluster", label=TRUE, repel = FALSE) + NoLegend()


out <- se.combined@meta.data
umap <- se.combined@reductions$umap@cell.embeddings
out <- cbind(out, umap)

#ggplot(out, aes(x=UMAP_1, y=UMAP_2, color=as.factor(queryCluster))) + geom_point(size=0.2) + theme_classic()
#ggplot(out[which(!is.na(out$refCluster)), ], aes(x=UMAP_1, y=UMAP_2, color=as.factor(refCluster))) + geom_point(size=0.2) + theme_classic()

#------------------------
# overlap score
metaout <- out

#-----------------------------------------------------
# t1: table with 2 columns: coembed labels, raw labels
# t2: table with 2 columns: coembed labels, raw labels
cal_ovlpScore <- function(t1, t2){
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- plyr::join(t1.pct.df, t2.pct.df, by="ident", type="full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}


# calculate overlap
ident2ref <- data.frame(idents=metaout$seurat_clusters, rna_label=metaout$refCluster)
ident2ref <- ident2ref[complete.cases(ident2ref), ]

ident2query <- data.frame(idents=metaout$seurat_clusters, atac_label=metaout$queryCluster)
ident2query <- ident2query[complete.cases(ident2query), ]

ovlpScore.df <- cal_ovlpScore(ident2ref, ident2query)
ovlpScore.df <- as.data.frame(ovlpScore.df)
#write.table(ovlpScore.df, "./rds/mba.hba.coembed.overlapScore.tsv", sep="\t", col.names = T, row.names = F, quote = F)

# ovlpScore.df.sel <- subset(ovlpScore.df, ovlpScore.df$ovlpScore>=0.2)
ovlpScore.mx <- reshape2::dcast(ovlpScore.df, anno1~anno2, value.var="ovlpScore", fill = 0)
ovlpScore.mx <- as.data.frame(ovlpScore.mx)
#write.table(ovlpScore.mx, "./rds/mba.hba.coembed.overlapScore.mx.tsv", sep="\t", col.names = T, row.names = F, quote = F)

ovlpScore.plot <- ovlpScore.mx
rownames(ovlpScore.plot) <- ovlpScore.plot$anno1
ovlpScore.plot$anno1 <- NULL
ovlpScore.plot <- as.matrix(ovlpScore.plot)

pdf(paste0(outF, ".combinedPlot.pdf"))
pt_combined_elbow
pt_coembed
pt_umap_refCluster
pt_umap_queryCluster
pt_umap_combineCluster
pheatmap(ovlpScore.plot, main = "overlap score: row: Allen; col: Aging", fontsize_row = 5, fontsize_col = 5)
dev.off()







