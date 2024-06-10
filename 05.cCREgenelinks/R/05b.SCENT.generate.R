#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
#parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-n", "--nbatch", required=TRUE, help="nbatch of peak list")
parser$add_argument("-d1", "--dir1", required=TRUE, help="directory to the project")
parser$add_argument("-d2", "--dir2", required=TRUE, help="directory to the project")
parser$add_argument("-s", "--subclass", required=TRUE, help="subclass to subset")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

#input = as.character(args$input)
nbatch = as.numeric(args$nbatch)
path1 = as.character(args$dir1)
path2 = as.character(args$dir2)
subclass = as.character(args$subclass)
outF = as.character(args$output)

suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
suppressMessages(library("SCENT"))
options(stringsAsFactors = FALSE)

all_gene_peak <- readRDS(file=file.path(path1, "rds/ATAC/after_integra/GRN/SCENT/data/SN_integra.gene_peak_pairs.rds"))
SCENT_obj <- readRDS(file=file.path(path1, "rds/ATAC/after_integra/GRN/SCENT/data/SN_integra.SCENT_obj.rds"))

split_SCENT_batch <- function(object, gene_peak, nbatch=1000) {
  gene_peak$index <- 1:nrow(gene_peak)
  gene_peak$batch_index <- Hmisc::cut2(gene_peak$index, g = nbatch, levels.mean = TRUE)
  gene_peak_list <- split(gene_peak, f = gene_peak$batch_index)
  gene_peak_list <- lapply(gene_peak_list, function(x) x[(names(x) %in% c("peak", "gene"))])
  names(gene_peak_list) <- 1:length(gene_peak_list)
  # Update the SCENT.peak.info field of the constructor in R:
  object@peak.info.list <- gene_peak_list
  return(object)
}

subclass_rename <- gsub(" ", "_", subclass)
subclass_peaks <- read.table(paste0(path2, "/peak_calling/after_integra/final/subclass.final.peak.srt/", subclass_rename, ".bed"))
subclass_called_peaks <- paste0(subclass_peaks$V1,":",subclass_peaks$V2,"-",subclass_peaks$V3)
message("!!!!!", subclass, " called peaks number: ", length(subclass_called_peaks), "!!!!!")
subclass_gene_peak <- all_gene_peak[which(all_gene_peak$peak %in% subclass_called_peaks),]
message("!!!!!", subclass, " gene-peak pairs number: ", nrow(subclass_gene_peak), "!!!!!")
subclass_SCENT_obj <- split_SCENT_batch(SCENT_obj, subclass_gene_peak, nbatch=nbatch) 

saveRDS(subclass_SCENT_obj, file=file.path(outF, paste0(subclass_rename, ".SCENT.nbatch", nbatch, ".rds")))

print("Job is done!")


