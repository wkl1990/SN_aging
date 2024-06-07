#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input pca rds")
parser$add_argument("-a", "--assay", default="integrated", help="assay")
parser$add_argument("-m", "--method", default="all", help="methods")
parser$add_argument("-c", "--npc", default=NULL, help="number of pcs to estimate")
parser$add_argument("-t", "--plot", action='store_true', help="whether to output plot")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
assay = as.character(args$assay)
if (!is.null(args$npc)){
  npc = as.character(args$npc)
}
method = as.character(args$method)
plot = ifelse(args$plot, TRUE, FALSE) 
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("findPC"))


autoPC_rna <- function(RNAdata, npc=100, assay="RNA", method="all") {
  DefaultAssay(RNAdata) <- assay
  data.use <- Stdev(RNAdata)
  pc_len <- length(data.use)
  npc <- ifelse(pc_len>=max(npc), npc, unique(ceiling(seq(pc_len, pc_len/3, length.out=length(npc)))))
  pc_mean <- findPC(sdev=data.use, number=npc, method=method, aggregate='mean')
  pc_median <- findPC(sdev=data.use, number=npc, method=method, aggregate='median')
  pc_voting <- findPC(sdev=data.use, number=npc, method=method, aggregate='voting')
  pc_all <- findPC(sdev=data.use, number=npc, method=method)
  pc_combine <- data.frame(mean=pc_mean, median=pc_median, voting=pc_voting)
  pclist <- list(data=data.use, npc=npc, all=pc_all, combine=pc_combine)
  return(pclist)
}

print("Step1: read rds data!")
RNAdata <- readRDS(input)
if (exists(quote(npc))) {
  npc <- as.numeric(unlist(str_split(npc, ",")))
} else {
  print("Number of pcs not provided, use seq(30,100,10) default npc!")
  npc <- seq(30, 100, 10)
}

print("Step2: pc estimate start!")
if (length(RNAdata)==1) {
  RNA_pc <- autoPC_rna(RNAdata, npc=npc, assay=assay, method=method)
  data.use <- RNA_pc$data
  pc_all <- RNA_pc$all
  pc_combine <- RNA_pc$combine 
  npc <- RNA_pc$npc
} else {
  pc_all <- list()
  pc_combine <- list()
  data.use <- list()
  npc_list <- list()
  for (i in names(RNAdata)) {
    RNA_pc <- autoPC_rna(RNAdata[[i]], npc=npc, assay=assay, method=method)
    pc_all[[i]] <- RNA_pc$all
    pc_combine[[i]] <- RNA_pc$combine
    data.use[[i]] <- RNA_pc$data
    npc_list[[i]] <- RNA_pc$npc
    rm(RNA_pc)      
  }
}

print("Step3: save data!")
if (length(RNAdata)==1) {
  if (plot==TRUE){
    print("Output findPC plot!")
    pdf(paste0(outF, ".findPC.pdf"), width=16, height=8)
    findPC(sdev=data.use, number=npc, method=method, figure=TRUE)
    dev.off()
  }
  write.csv(pc_all, file=paste0(outF, ".pc_all.csv"))
  write.csv(pc_combine, file=paste0(outF, ".pc_combine.csv"))
} else {
  if (plot==TRUE){
    print("Output findPC plot!")
    pdf(paste0(outF, ".findPC.pdf"), width=16, height=8)
    for (i in names(RNAdata)) {
      findPC(sdev=data.use[[i]], number=npc_list[[i]], method=method, figure=TRUE)
    }
    dev.off()
  }
  for (i in names(RNAdata)) {
    write.csv(pc_all, file=paste0(outF, ".", i, ".pc_all.csv"))
    write.csv(pc_combine, file=paste0(outF, ".", i, ".pc_combine.csv"))
  }
}
print("Job is done!")



