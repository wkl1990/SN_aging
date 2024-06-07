#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input h5 file(s)")
parser$add_argument("-d", "--path", default="./", help="directory path")
parser$add_argument("-n", "--name", default=NULL, help="sample name")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
path = as.character(args$path)
if (!is.null(args$name)) {
  name = as.character(args$name)
}
outF = as.character(args$output)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("stringr"))

# Load one dataset
load_1data <- function(folder, sample=NULL){
  if (is.null(sample)) {
    sample <- folder
  }
  names(folder) <- sample
  print(paste0("Load ", sample, " start."))
  data <- Read10X_h5(folder)
  #extract RNA data
  data_rna_counts <- data$`Gene Expression`
  #Create Seurat object
  RNAdata <- CreateSeuratObject(counts=data_rna_counts, project=sample, min.cells=0, min.features=0)
  RNA_cells <- paste(sample, colnames(RNAdata), sep=":")
  RNAlist <- list(data=RNAdata, cell=RNA_cells)
  return(RNAlist)
}


# Load multiple samples
load_2data <- function(folders, samples=NULL){
  if (is.null(samples)) {
    samples <- folders
  }
  names(folders) <- samples
  RNAsceList <- list()
  for (sample in samples){
    print(paste0("Load ", sample, " start."))
    RNAsceList[[sample]] <- CreateSeuratObject(counts=Read10X_h5(folders[sample])$`Gene Expression`, project=sample, min.cells=0, min.features=0)
  } 
  RNA_cells_raw <- c()
  for (i in names(RNAsceList)){
    RNA_cells_raw <- c(RNA_cells_raw, paste(i, colnames(RNAsceList[[i]]), sep=":"))
  }
  RNAlist <- list(data=RNAsceList, cell=RNA_cells_raw)
  return(RNAlist)
}

print("Step1: parse argument!")
paths <- unlist(str_split(path, ","))
files <- unlist(str_split(input, ","))
folders <- paste0(paths, files)

print("Step2: load data!")
if (exists(quote(name))) {
  samples <- unlist(str_split(name, ","))
  if (length(files)!=length(samples)) {
    stop("Length of sample names does not match files!")    
  } else {
    if (length(folders)==1) {
      RNAlist <- load_1data(folders, samples)
    } else {
      RNAlist <- load_2data(folders, samples)
    } 
  }
} else {
  print("No sample names are provided, will use file names!")
  if (length(folders)==1) {
    RNAlist <- load_1data(folders)
  } else {
    RNAlist <- load_2data(folders)
  } 
}
RNAdata <- RNAlist$data
RNAcell <- RNAlist$cell

print("Step3: save data!")
saveRDS(RNAdata, file=paste0("./rds/", outF, ".raw.rds"))
saveRDS(RNAcell, file=paste0("./rds/", outF, ".cells_raw.rds"))
print("Job is done!")

