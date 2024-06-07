#!/usr/bin/env Rscript

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="input pca rds")
parser$add_argument("-a", "--assay", default="integrated", help="assay")
parser$add_argument("-m", "--matrix", default="pca", help="which matrix to calculate silhouette")
#parser$add_argument("-c", "--cluster", default="seurat_clusters", help="cluster info in the meta data")
#parser$add_argument("-r", "--recluster", action="store_true", help="do the clustering")
parser$add_argument("-t", "--plot", action='store_true', help="whether to output plot")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

input = as.character(args$input)
assay = as.character(args$assay)
matrix = as.character(args$matrix)
#cluster = as.character(args$cluster)
#recluster = ifelse(args$recluster, TRUE, FALSE) 
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
suppressPackageStartupMessages(library("cluster"))
suppressPackageStartupMessages(library("parallelDist"))


#cal_silhouette <- function(x, cluster){
#  if(unique(cluster) < 2) stop("Cluster level must bet > = 2")
#  if (is.data.frame(x)) x <- as.matrix(x)
#  if(is.null(diss)) diss <- stats::dist(x)
#  ss <- cluster::silhouette(cluster, diss)
#  return(ss)  
#}

cal_silhouette <- function(x, cluster){
  if(length(unique(cluster)) < 2) {
    message("Only one cluster level were found!")
    return(NA)
  }
  diss <- parDist(x)
  sil <- cluster::silhouette(cluster, diss)
  ss <- mean(sil[,3]) 
  return(ss)  
}

cluster_rna <- function(RNAdata, resolution=resolution, algorithm=method, assay="integrated") {
  DefaultAssay(RNAdata) <- assay
  RNAdata <- FindClusters(object=RNAdata, resolution=resolution, algorithm=algorithm)
  return(RNAdata)
}

print("Step1: read rds data!")
RNAdata <- readRDS(input)
DefaultAssay(RNAdata) <- assay

if (matrix=="pca") {
  RNA_mat <- Embeddings(RNAdata, reduction="pca")
} else {
  RNA_data <- GetAssayData(RNAdata, slot="scale.data")
  RNA_mat <- t(RNA_data)
}

print("Step2: run clustering and silhouette!")
umap_list <- list()
resolutions <- seq(0.1,2,0.1)
silhouette_tbl <- data.frame(Resolution=resolutions, Sil_score=rep(NA,length(resolutions)))
for (resolution in resolutions) {
  message("Calculate silhouette score for resolution ", resolution)
  RNAdata <- cluster_rna(RNAdata, resolution=resolution, algorithm=4, assay=assay)
  RNA_cluster <- as.integer(RNAdata$seurat_clusters)
  sil_score <- cal_silhouette(RNA_mat, RNA_cluster)
  umap_list[[as.character(resolution)]] <- DimPlot(RNAdata, label=TRUE, repel=FALSE) + NoLegend() + ggtitle(paste0("resolution: ", resolution))
  silhouette_tbl$Sil_score[which(silhouette_tbl$Resolution==resolution)] <- sil_score
}
silhouette_select <- silhouette_tbl$Resolution[which.max(silhouette_tbl$Sil_score)]

print("Step3: output silhouette score!")
write.csv(silhouette_tbl, file=paste0(outF, ".silhouette.csv"), row.names=FALSE, quote=FALSE)
write.table(silhouette_select, file=paste0(outF, ".resolution_select.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
if (plot) {
  pdf(file = paste0(outF, ".silhouette.pdf"), width = 15, height = 12)
  plot(ggarrange(plotlist=umap_list, nrow=5,ncol=4))
  plot(silhouette_tbl$Resolution, silhouette_tbl$Sil_score,
       type = "l", pch = 16, col = "red", lwd = 4,
       xlab = "Resolution", ylab = "Silhouette", cex.lab = 1.5,
       main = paste("Silhouette analysis for Leiden-base clustering")
       )
  dev.off()
}




