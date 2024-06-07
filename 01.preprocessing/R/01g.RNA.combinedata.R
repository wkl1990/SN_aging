#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input RNA rds")
parser$add_argument("-n", "--nfeature", default=2000, help="number of selected feature")
parser$add_argument("-c", "--npc", default=NULL, help="number of pc used")
parser$add_argument("-s", "--samples", default=NULL, help="samples merge info that joint by @")
parser$add_argument("-m", "--mincells", default=20, help="merge samples if samples contain less than the cutoff cells")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
nfeature = as.numeric(args$nfeature)
if (!is.null(args$npc)){
  npc = as.numeric(args$npc)
}
if (!is.null(args$samples)){
  samples = as.character(args$samples)
}
mincells = as.numeric(args$mincells)
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("stringr"))

print("Step1: read rds data!")
RNAdata <- readRDS(input)
if (exists(quote(npc))) {
  npc <- as.numeric(npc)
} else {
  print("Number of pcs not provided, use 20 default npc!")
  npc <- 20
}
if (exists(quote(samples))) {
  samples <- as.character(unlist(str_split(samples, ",")))
  sample_list <- str_split(samples, "@")
} else {
  print("samples merge not provided, be careful for small samples size!")
}

print("Step2: find anchors and integrate data!")
#normalize and identify variable features for each dataset independently
min_sample <- min(unlist(lapply(RNAdata,ncol)))
while (min_sample==1) {
  print(paste0("Only 1 cell in sample ", names(RNAdata)[which(unlist(lapply(RNAdata,ncol))==1)], "!"))
  RNAdata_merge <- list()
  for (i in 1:length(sample_list)) {
    if (length(sample_list[[i]])==2) {
      print(paste0("Too little cells in samples, merge samples for ", paste(sample_list[[i]], collapse=" "), "!"))
      RNAdata_merge[[i]] <- merge(RNAdata[[sample_list[[i]][1]]], RNAdata[[sample_list[[i]][2]]]) #, add.cell.ids=c("2m_rep1", "2m_rep2")
    } else if (length(sample_list[[i]])==1) {
      print(paste0("No replication to merge for sample ", paste(sample_list[[i]], collapse=" "), "!"))
      RNAdata_merge[[i]] <- RNAdata[[sample_list[[i]][1]]]
    } else {
      print("Too little cells in samples, try to merge samples but only support merge two samples now!")
    }
  }
  RNAdata <- RNAdata_merge
  min_sample <- min(unlist(lapply(RNAdata,ncol)))
  sample_list <- split(1:length(samples), ceiling(seq_along(1:length(samples))/2))
  samples <- 1:length(sample_list)
}

RNAdata <- lapply(X=RNAdata, FUN=function(RNAsce) {
    RNAsce <- NormalizeData(RNAsce, verbose=FALSE)
    RNAsce <- FindVariableFeatures(RNAsce, selection.method="vst", nfeatures=nfeature, verbose=FALSE)
})


#select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list=RNAdata, nfeatures=nfeature, verbose=FALSE)

#Integration
min_sample <- min(unlist(lapply(RNAdata,ncol)))
while (min_sample<mincells) {
  RNAdata_merge <- list()
  for (i in 1:length(sample_list)) {
    if (length(sample_list[[i]])==2) {
      print(paste0("Too little cells in samples, merge samples for ", paste(sample_list[[i]], collapse=" "), "!"))
      RNAdata_merge[[i]] <- merge(RNAdata[[sample_list[[i]][1]]], RNAdata[[sample_list[[i]][2]]]) #, add.cell.ids=c("2m_rep1", "2m_rep2")
    } else if (length(sample_list[[i]])==1) {
      print(paste0("No replication to merge for sample ", paste(sample_list[[i]], collapse=" "), "!"))
      RNAdata_merge[[i]] <- RNAdata[[sample_list[[i]][1]]]
    } else {
      print("Too little cells in samples, try to merge samples but only support merge two samples now!")
    }
  }
  RNAdata <- RNAdata_merge
  sample_list <- split(1:length(samples), ceiling(seq_along(1:length(samples))/2))
  samples <- 1:length(sample_list)
  min_sample <- min(unlist(lapply(RNAdata,ncol)))
  if (length(RNAdata)<2) {
    break
  }
}   
#npc <- ifelse(min_sample>npc*2, npc, min_sample/2)
if (length(RNAdata)<2) {
  print("Too little cells in samples, merge all samples instead of integration!")
  print("Step3: save data!")
  write.csv(lapply(RNAdata,dim), file=paste0(outF, ".merge_samples_table.csv"))
  saveRDS(RNAdata[[1]], file=paste0(outF, ".merged.rds"))
} else {
  RNA.anchors <- FindIntegrationAnchors(object.list=RNAdata, anchor.features=features, dims=1:npc, verbose=FALSE)
  k.weight <- ifelse(min_sample>100 & min_sample>mincells*2, 100, floor(min_sample/2)-1)
  RNA.combined <- IntegrateData(anchorset=RNA.anchors, dims=1:npc, k.weight=k.weight, verbose=FALSE)

  print("Step3: save data!")
  write.csv(lapply(RNAdata,dim), file=paste0(outF, ".merge_samples_table.csv"))
  saveRDS(RNA.combined, file=paste0(outF, ".combined.rds"))
}

print("Job is done!")


