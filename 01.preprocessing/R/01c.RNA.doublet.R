#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input cluster rds")
parser$add_argument("-c", "--npc", default=NULL, help="number of pc used")
parser$add_argument("-m", "--mincells", default=1, help="genes detected in at least the cutoff cells")
parser$add_argument("-f", "--filter", action='store_true', help="whether to filter doublet")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
if (!is.null(args$npc)){
  npc = as.character(args$npc)
}
mincells = as.numeric(args$mincells)
dfilter = ifelse(args$filter, TRUE, FALSE) 
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("DoubletFinder"))

doublet_rna <- function(RNAdata, npc=20, filter=TRUE, min.cells=1) {
  doublet_rate <- data.frame(
    rate=c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
    loadedCells=c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
    recoverCells=c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
  # pre-process (pca-umap-cluster) should be finished for RNAdata
  # pK Identification
  set.seed(2023)
  sweep.res.list_RNAdata <- paramSweep_v3(RNAdata, PCs=1:npc, sct=FALSE)
  sweep.stats_RNAdata <- summarizeSweep(sweep.res.list_RNAdata, GT=FALSE)
  bcmvn_RNAdata <- find.pK(sweep.stats_RNAdata)
  # Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- RNAdata@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- RNAsceList[[1]]@meta.data$ClusteringResults
  recoverNum <- round(ncol(RNAdata)/1000)*1000
  if(recoverNum<500){recoverNum=500}
  if(recoverNum>10000){recoverNum=10000}
  print(paste0("Estimated cell number is: ", recoverNum))
  dourate <- doublet_rate[match(recoverNum, doublet_rate$recoverCells), "rate"] ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi <- round(dourate*nrow(RNAdata@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  # Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  bcmvn_RNAdata <- bcmvn_RNAdata[order(bcmvn_RNAdata$BCmetric,decreasing=TRUE),]
  pK <- as.numeric(as.character(bcmvn_RNAdata[1,"pK"]))
  print(paste0("Estimated optimal pK is: ", pK))
  RNAdata <- doubletFinder_v3(RNAdata, PCs=1:npc, pN=0.25, pK=pK, nExp=nExp_poi.adj, reuse.pANN=FALSE, sct=FALSE)
  # filter
  if (filter==TRUE) {
    RNAdata <- RNAdata[, which(RNAdata@meta.data[ ,grep("DF", colnames(RNAdata@meta.data))] == "Singlet")]
    RNAdata_counts <- GetAssayData(RNAdata, slot="count")
    RNAdata_genes <- which(rowSums(RNAdata_counts>0)>=min.cells)
    #print(length(RNAdata_genes))
    RNAdata <- RNAdata[RNAdata_genes,]    
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

print("Step2: DoubletFinder start!")
if (length(RNAdata)==1) {
  RNAdata_doubletlist <- doublet_rna(RNAdata, npc=npc, filter=dfilter, min.cells=mincells)
  RNAsceList <- RNAdata_doubletlist
  RNA_cells_doublet <- paste(names(RNAdata), colnames(RNAdata), sep=":")
} else {
  RNAsceList <- list()
  RNA_cells_doublet <- c()
  for (i in names(RNAdata)) {
    RNAdata_doubletlist <- doublet_rna(RNAdata[[i]], npc=npc[[i]], filter=dfilter, min.cells=mincells)
    RNAsceList[[i]] <- RNAdata_doubletlist
    RNA_cells_doublet <- c(RNA_cells_doublet, paste(i, colnames(RNAsceList[[i]]), sep=":"))
    rm(RNAdata_doubletlist)
  }
}

print("Step3: save data!")
saveRDS(RNAsceList, file=paste0("./rds/", outF, ".doublet.rds"))
saveRDS(RNA_cells_doublet, file=paste0("./rds/", outF, ".cells_doublet.rds"))
print("Job is done!")


