#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
#parser$add_argument("-i", "--input", required=TRUE, help="input RDS")
parser$add_argument("-n", "--node", required=TRUE, help="node of peak list")
parser$add_argument("-c", "--core", default=4, help="cores")
parser$add_argument("-d", "--dir", required=TRUE, help="directory to the project")
parser$add_argument("-s", "--subclass", required=TRUE, help="subclass to subset")
parser$add_argument("-r", "--regr", default="poisson", help="Regression Type")
#parser$add_argument("-b", "--bin", action='store_true', help="Binarize ATAC counts")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

#input = as.character(args$input)
node = as.numeric(args$node)
cores = as.numeric(args$core)
path = as.character(args$dir)
subclass = as.character(args$subclass)
regr = as.character(args$regr)
#bin = ifelse(args$plot, TRUE, FALSE) 
bin = TRUE
outF = as.character(args$output)

suppressMessages(library("SCENT"))
options(stringsAsFactors = FALSE)

subclass_rename <- gsub(" ", "_", subclass)
SCENT_obj <- readRDS(file=file.path(path, "rds/ATAC/after_integra/GRN/SCENT/data/", paste0(subclass_rename, ".SCENT.nbatch1000.rds")))

#### Get the corresponding dataframe from the list:
SCENT_obj@peak.info <- SCENT_obj@peak.info.list[[node]]

#### Run SCENT algorithm of Tnk cell type and use 6 cores for parallelization:
SCENT_obj <- SCENT_algorithm(SCENT_obj, subclass, cores, regr, bin)

#### Output SCENT results for each gene-peak pair block.
filename <- file.path(outF, paste0(subclass_rename, ".SCENTresult.node", node, ".txt")) 

if (!dir.exists(outF)) {
  dir.create(outF, recursive = TRUE)
}
write.table(SCENT_obj@SCENT.result, file = filename, row.names = F, col.names = T)

print("Job is done!")


