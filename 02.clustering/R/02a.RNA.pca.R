#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input qc rds")
parser$add_argument("-a", "--assay", default="integrated", help="assay")
parser$add_argument("-n", "--nFeature", default=2000, help="number of variable features")
parser$add_argument("-c", "--npc", default=100, help="number of pcs to estimate")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
assay = as.character(args$assay)
nFeature = as.numeric(args$nFeature)
npc = as.numeric(args$npc)
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))

preprocess_rna <- function(RNAdata, nFeature=2000, npc=100, assay="RNA") {
  DefaultAssay(RNAdata) <- assay
  RNAdata <- NormalizeData(RNAdata, verbose=FALSE)
  RNAdata <- FindVariableFeatures(RNAdata, selection.method="vst", nfeatures=nFeature, verbose=FALSE)
  RNAdata <- ScaleData(RNAdata, verbose=FALSE)
  cellnum <- ncol(RNAdata)
  npc <- ifelse(cellnum>npc, npc, cellnum-1)
  RNAdata <- RunPCA(RNAdata, npcs=npc, verbose=FALSE)
  pt_elbow <- ElbowPlot(RNAdata, ndims=npc)
  pclist <- list(data=RNAdata, pcplot=pt_elbow)
  return(pclist)
}

preprocess_integrated_rna <- function(RNAdata, npc=100, assay="integrated") {
  DefaultAssay(RNAdata) <- assay
  RNAdata <- ScaleData(RNAdata, verbose=FALSE)
  cellnum <- ncol(RNAdata)
  npc <- ifelse(cellnum>npc, npc, cellnum-1)
  RNAdata <- RunPCA(RNAdata, npcs=npc, verbose=FALSE)
  pt_elbow <- ElbowPlot(RNAdata, ndims=npc)
  pclist <- list(data=RNAdata, pcplot=pt_elbow)
  return(pclist)
}

print("Step1: read rds data!")
RNAdata <- readRDS(input)

print("Step2: pc estimate start!")
if (length(RNAdata)==1) {
  if (assay=="RNA") {
    RNAdata_pclist <- preprocess_rna(RNAdata, nFeature=nFeature, npc=npc)
    RNAsceList <- RNAdata_pclist$data
    pt_elbow <- RNAdata_pclist$pcplot    
  } else if (assay=="integrated") {
    RNAdata_pclist <- preprocess_integrated_rna(RNAdata, npc=npc)
    RNAsceList <- RNAdata_pclist$data
    pt_elbow <- RNAdata_pclist$pcplot
  }
} else {
  pt_elbow <- list()
  RNAsceList <- list()
  for (i in names(RNAdata)) {
    if (assay=="RNA") {
      RNAdata_pclist <- preprocess_rna(RNAdata[[i]], nFeature=nFeature, npc=npc)
      RNAsceList[[i]] <- RNAdata_pclist$data
      pt_elbow[[i]] <- RNAdata_pclist$pcplot
      rm(RNAdata_pclist)      
    } else if (assay=="integrated") {
      RNAdata_pclist <- preprocess_integrated_rna(RNAdata[[i]], npc=npc)
      RNAsceList[[i]] <- RNAdata_pclist$data
      pt_elbow[[i]] <- RNAdata_pclist$pcplot
      rm(RNAdata_pclist)      
    }
  }
}

print("Step3: save data!")
if (length(RNAdata)==1) {
  ggsave(pt_elbow, file=paste0(outF, ".pc_elbow.pdf"), width=16, height=8)
} else {
  for (i in names(pt_elbow)) {
    ggsave(pt_elbow[[i]], file=paste0(outF, ".", i, ".pc_elbow.pdf"), width=16, height=8)
  }
}
saveRDS(RNAsceList, file=paste0(outF, ".pca.rds"))
print("Job is done!")


