#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input raw rds")
parser$add_argument("-n", "--nFeature", default=500, help="nFeature cutoff")
parser$add_argument("-m", "--mt", default=5, help="mitochondria cutoff")
parser$add_argument("-c", "--mincells", default=1, help="genes detected in at least the cutoff cells")
parser$add_argument("-g", "--organism", default="mouse", help="organism species")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
organism = as.character(args$organism)
nFeature = as.numeric(args$nFeature)
mt = as.numeric(args$mt)
mincells = as.numeric(args$mincells)
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))


qc_rna <- function(RNAdata, organism="mouse", nFeature=500, mt=5, min.cells=1) {
  if (organism=="mouse") {
    pattern <- "^mt-"
  } else if (organism=="human") {
    pattern <- "^MT-"
  } else {
    stop("Unrecognized organism!")
  }
  RNAdata[["percent.mt"]] <- PercentageFeatureSet(RNAdata, pattern=pattern)
  RNAdata <- subset(RNAdata, subset=nFeature_RNA>nFeature & percent.mt<mt)
  RNAdata_counts <- GetAssayData(RNAdata, slot="count")
  RNAdata_genes <- which(rowSums(RNAdata_counts>0)>=min.cells)
  #print(length(RNAdata_genes))
  RNAdata <- RNAdata[RNAdata_genes,]    
  pt_qc_vln <- VlnPlot(RNAdata, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
  pt_qc_UMIvsMT <- FeatureScatter(RNAdata, feature1="nFeature_RNA", feature2="percent.mt")
  qclist <- list(data=RNAdata, vlnplot=pt_qc_vln, scatplot=pt_qc_UMIvsMT)
  return(qclist)
}

print("Step1: read rds data!")
RNAdata <- readRDS(input)

# QC
print("Step2: qc start!")
if (length(RNAdata)==1) {
  RNAdata_qclist <- qc_rna(RNAdata, organism=organism, nFeature=nFeature, mt=mt, min.cells=mincells)
  RNAsceList <- RNAdata_qclist$data
  pt_qc_vln <- RNAdata_qclist$vlnplot
  pt_qc_UMIvsMT <- RNAdata_qclist$scatplot
  RNA_cells_qc <- paste(names(RNAdata), colnames(RNAdata), sep=":")
} else {
  pt_qc_vln <- list()
  pt_qc_UMIvsMT <- list()
  RNAsceList <- list()
  RNA_cells_qc <- c()
  for (i in names(RNAdata)) {
    RNAdata_qclist <- qc_rna(RNAdata[[i]], organism=organism, nFeature=nFeature, mt=mt, min.cells=mincells)
    RNAsceList[[i]] <- RNAdata_qclist$data
    pt_qc_vln[[i]] <- RNAdata_qclist$vlnplot
    pt_qc_UMIvsMT[[i]] <- RNAdata_qclist$scatplot
    RNA_cells_qc <- c(RNA_cells_qc, paste(i, colnames(RNAsceList[[i]]), sep=":"))
    rm(RNAdata_qclist)
  }
}

print("Step3: save data!")
if (length(RNAdata)==1) {
  ggsave(pt_qc_vln, file=paste0("./plot/", outF, ".qc.vln.pdf"), width=16, height=8)
  ggsave(pt_qc_UMIvsMT, file=paste0("./plot/", outF, ".qc.UMIvsMT.pdf"), width=8, height=8)
} else {
  for (i in names(pt_qc_vln)) {
    ggsave(pt_qc_vln[[i]], file=paste0("./plot/", outF, ".", i, ".qc.vln.pdf"), width=16, height=8)
    ggsave(pt_qc_UMIvsMT[[i]], file=paste0("./plot/", outF, ".", i, ".qc.UMIvsMT.pdf"), width=8, height=8)
  }
}
saveRDS(RNAsceList, file=paste0("./rds/", outF, ".qc.rds"))
saveRDS(RNA_cells_qc, file=paste0("./rds/", outF, ".cells_qc.rds"))
print("Job is done!")




