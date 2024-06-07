#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input umap RDS")
parser$add_argument("-r", "--reduction", default="pca", help="reduction to use")
parser$add_argument("-a", "--assay", default="integrated", help="assay")
parser$add_argument("-c", "--npc", default=NULL, help="number of pc used")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
reduction = as.character(args$reduction)
assay = as.character(args$assay)
if (!is.null(args$npc)){
  npc = as.character(args$npc)
}
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))

knn_rna <- function(RNAdata, reduction="pca", npc=20, assay="integrated") {
  DefaultAssay(RNAdata) <- assay
  RNAdata <- FindNeighbors(RNAdata, reduction=reduction, dims=1:npc, verbose=FALSE)
  return(RNAdata)
}

print("Step1: read rds data!")
RNAdata <- readRDS(input)
if (exists(quote(npc))) {
  npc <- as.numeric(unlist(str_split(npc, ",")))
} else {
  print("Number of pcs not provided, use 20 default npc!")
  npc <- rep(20, length(RNAdata))
}
if (length(npc)>1) {
  names(npc) <- names(RNAdata)
}

print("Step2: run KNN!")
if (length(RNAdata)==1) {
  RNAdata_knnlist <- knn_rna(RNAdata, reduction=reduction, npc=npc, assay=assay)
  RNAsceList <- RNAdata_knnlist
} else {
  RNAsceList <- list()
  for (i in names(RNAdata)) {
    RNAdata_knnlist <- knn_rna(RNAdata[[i]], reduction=reduction, npc=npc[[i]], assay=assay)
    RNAsceList[[i]] <- RNAdata_knnlist
    rm(RNAdata_knnlist)
  }
}


print("Step3: save data!")
knn_name <- paste0(assay, "_", "nn")
if (length(RNAdata)==1) {
  Matrix::writeMM(RNAsceList@graphs[[knn_name]], file=paste0(outF, ".knn.mmtx"))
} else {
  for (i in names(RNAdata)) {
    Matrix::writeMM(RNAsceList[[i]]@graphs[[knn_name]], file=paste0(outF, ".", i, ".knn.mmtx"))
  }
}

saveRDS(RNAsceList, file=paste0(outF, ".knn.rds"))
print("Job is done!")


