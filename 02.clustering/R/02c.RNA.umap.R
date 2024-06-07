#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input cluster rds")
parser$add_argument("-a", "--assay", default="integrated", help="assay")
parser$add_argument("-c", "--npc", default=NULL, help="number of pc used")
parser$add_argument("-u", "--umap", default="raw", help="set a and b to use raw umap")
parser$add_argument("-t", "--plot", action='store_true', help="whether to output plot")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
if (!is.null(args$npc)){
  npc = as.character(args$npc)
}
assay = as.character(args$assay)
umap = as.character(args$umap)
uplot = ifelse(args$plot, TRUE, FALSE) 
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))

umap_rna <- function(RNAdata, npc=20, umap="raw", assay="integrated") {
  DefaultAssay(RNAdata) <- assay
  if (umap=="raw") {
    if (ncol(RNAdata)<=30) {
      print("Warning: Sample number is less than 30!")
      RNAdata <- RunUMAP(RNAdata, dims=1:npc, a=1.8956, b=0.8006, verbose=FALSE, n.neighbors=ceiling(ncol(RNAdata)/2))
    } else {
      RNAdata <- RunUMAP(RNAdata, dims=1:npc, a=1.8956, b=0.8006, verbose=FALSE)
    }
  } else if (umap=="seurat") {
    if (ncol(RNAdata)<=30) {
      print("Warning: Sample number is less than 30!")
      RNAdata <- RunUMAP(RNAdata, dims=1:npc, verbose=FALSE, n.neighbors=ceiling(ncol(RNAdata)/2))
    } else {
      RNAdata <- RunUMAP(RNAdata, dims=1:npc, verbose=FALSE)
    }
  } else {
    stop("Unrecognized umap!")
  }
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


print("Step2: run umap start!")
if (length(RNAdata)==1) {
  RNAdata_umaplist <- umap_rna(RNAdata, npc=npc, umap=umap, assay=assay)
  RNAsceList <- RNAdata_umaplist
} else {
  RNAsceList <- list()
  for (i in names(RNAdata)) {
    RNAdata_umaplist <- umap_rna(RNAdata[[i]], npc=npc[[i]], umap=umap, assay=assay)
    RNAsceList[[i]] <- RNAdata_umaplist
    rm(RNAdata_umaplist)
  }
}

print("Step3: save data!")
if (uplot==TRUE){
  print("Output umap plot!")
  #pt_umap <- DimPlot(RNAsceList, reduction="umap")
  #ggsave(paste0("./plot/", outF, ".raw.umap.pdf"), plot=pt_umap, width=8, height=8)
}
saveRDS(RNAsceList, file=paste0(outF, ".umap.rds"))
print("Job is done!")


