#!/usr/bin/R

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-r", "--ref", required=TRUE, help="ref input RNA rds")
parser$add_argument("-q", "--query", required=TRUE, help="query input RNA rds")
parser$add_argument("-a1", "--ref_assay", required=TRUE, help="ref assay")
parser$add_argument("-a2", "--query_assay", required=TRUE, help="query assay")
parser$add_argument("-c1", "--ref_cluster", required=TRUE, help="ref cluster in meta table")
parser$add_argument("-c2", "--query_cluster", required=TRUE, help="query cluster in meta table")
#parser$add_argument("-n", "--normalized", default="Strandard", help="normalization method, Strandard or SCTransform")
#parser$add_argument("-i", "--integrated", default="cca", help="integrated method, cca, rpca, jpca, or rlsi")
parser$add_argument("-t", "--transfered", default="cca", help="transfered method, cca, rpca, or pcaproject")
parser$add_argument("-d", "--dynamic_downsample", action="store_true", help="dynamic downsampling number of cells for data")
parser$add_argument("-d1", "--ref_downsample", default=NULL, help="downsampling number of cells for ref data")
parser$add_argument("-d2", "--query_downsample", default=NULL, help="downsampling number of cells for query data")
parser$add_argument("-m", "--remove", default=NULL, help="filter celltype with small cell number in reference")
parser$add_argument("-k", "--kanchor", default=5, help="number of anchors")
#parser$add_argument("-f", "--nfeature", default=2000, help="number of selected feature")
parser$add_argument("-c", "--npc", default=25, help="number of pc used")
parser$add_argument("-s", "--score", action="store_true", default=FALSE, help="output consensus score and unasigned clusters")
parser$add_argument("-f", "--cutoff", default=NULL, help="cutoff of unasigned clusters")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

ref = as.character(args$ref)
query = as.character(args$query)
ref_assay = as.character(args$ref_assay)
query_assay = as.character(args$query_assay)
ref_cluster = as.character(args$ref_cluster)
query_cluster = as.character(args$query_cluster)
#normalized = as.character(args$normalized)
transfered = as.character(args$transfered)
wet_reduc <- ifelse(transfered == "cca", "cca", "pca")
if(transfered == "cca") {
    ref_reduc <- NULL} else {
        ref_reduc <- "pca"}
dynamic_downsample = ifelse(args$dynamic_downsample, TRUE, FALSE)
if(!is.null(args$ref_downsample)) {ref_downsample = as.numeric(args$ref_downsample)}
if(!is.null(args$query_downsample)) {query_downsample = as.numeric(args$query_downsample)}
if(!is.null(args$remove)) {remove_num = as.numeric(args$remove)}
kanchor = as.numeric(args$kanchor)
#nfeature = as.numeric(args$nfeature)
npc = as.numeric(args$npc)
out_score = ifelse(args$score, TRUE, FALSE) 
cutoff = ifelse(!is.null(args$cutoff), as.numeric(args$cutoff), 0.8)
outF = as.character(args$output)

#suppressPackageStartupMessages(library("SnapATAC"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("circlize"))
library("dendsort")
library("cowplot")
library("patchwork")


# * load rds 
ref_data <- readRDS(ref)
query_data <- readRDS(query)
ref_data@active.assay <- ref_assay
query_data@active.assay <- query_assay

SetAssayNull <- function(data) {
  assays <- setdiff(names(data@assays), data@active.assay)
  if (length(assays) != 0) {
    for (assay in assays) {
      data[[assay]] <- NULL
    }
  }
  return(data)
} 
ref_data <- SetAssayNull(ref_data)
query_data <- SetAssayNull(query_data)

if(!(ref_cluster %in% names(ref_data@meta.data))) {
  ref_data@meta.data[,ref_cluster] <- Idents(ref_data)
}
if(!(query_cluster %in% names(query_data@meta.data))) {
  query_data@meta.data[,query_cluster] <- Idents(query_data)
}


removeSeurat <- function(seu, groupBy = "cluster_id", remove = 10) {
  groups <- unlist(seu[[groupBy]])
  groups_table <- table(groups)
  groups_rm <- names(groups_table)[groups_table<remove]
  cells <- colnames(seu)[!(groups %in% groups_rm)]
  #seu_rm <- subset(x=seu, subset=groupBy %in% groups_rm, invert=TRUE)
  seu_rm <- subset(x=seu, cells=cells)
  return(seu_rm)
}

if (exists("remove_num")) {
  ref_data <- removeSeurat(ref_data, groupBy = ref_cluster, remove = remove_num)
#  query_data <- removeSeurat(query_data, groupBy = query_cluster, remove = remove_num)
}


get.downsample.fun <- function(minNum = 500, maxNum = 1000) {
  f <- function(index, labels) {
    set.seed(2024)
    stat <- table(labels)
    ## put divided in the first since sometimes the full production may
    ## introduce overflow problem in R since
    # only 32-bit integer are supported
    idealNum <- (stat / sum(stat)) * minNum * length(stat)
    nms <- names(idealNum)
    # FIXME: logic is too complex
    lindex <- lapply(nms, function(i) {
      curIndex <- index[which(labels %in% i)]
      if ( (stat[i] <= max(minNum, idealNum[i])) ) {
        return(curIndex)
      }
      tmpIndex <- sample(curIndex, size = idealNum[i], replace = FALSE)
      tmp2Index <- tmpIndex
      if(length(tmpIndex) <= minNum) {
        tmp2Index <- sample(curIndex, size = minNum, replace = FALSE)
      }
      if (length(tmpIndex) > maxNum) {
        tmp2Index <- sample(tmpIndex, size = maxNum, replace = FALSE)
      }
      return(tmp2Index)
    })
    return(sort(unlist(lindex)))
  }
  return(f)
}

downsampleSeurat <- function(seu,
                             groupBy = "cluster_id",
                             minNum = 1000,
                             maxNum = 2000) {
  fn.dp <- get.downsample.fun(
    minNum = minNum, maxNum = maxNum)
  allcells <- colnames(seu)
  ## NOTE: seu[[groupBy]] will return data.frame with one column
  labels <- seu@meta.data[[groupBy]]
  if(is.null(labels)) {
    warning(groupBy, " is not int meta.data.")
    return(seu)
  }
  dp.cells <- fn.dp(index = allcells, labels = labels)
  subset(seu, cells = dp.cells)
}

if (dynamic_downsample) {
  ref_cellnum <- ncol(ref_data)
  query_cellnum <- ncol(query_data)
  ref_stat <- table(ref_data@meta.data[ref_cluster])
  query_stat <- table(query_data@meta.data[query_cluster])
  if (ref_cellnum<=query_cellnum) {
    ref_downsample <- round(median(ref_stat), 0)
    query_downsample <- round(ref_downsample*length(ref_stat)/length(query_stat),0)
  } else {
    query_downsample <- round(median(query_stat), 0)
    ref_downsample <- round(query_downsample*length(query_stat)/length(ref_stat),0)
  }
  message("Dynamic downsample: ", "ref_downsample: ", ref_downsample, ", query_downsample: ", query_downsample)
  message("Number of clusters: ", "ref clusters: ", length(unique(ref_data@meta.data[[ref_cluster]])), ", query clusters: ", length(unique(query_data@meta.data[[query_cluster]])))  
}

if (exists("ref_downsample")) {
  #message("require: downsample Seurat.")
  #message(str_glue("atac: {atac_dsMin} - {atac_dsMax}."))
  ref_data <- downsampleSeurat(seu = ref_data, groupBy = ref_cluster,
    minNum = ref_downsample, maxNum = ref_downsample)
}
if (exists("query_downsample")) {
  #message("require: downsample Seurat.")
  #message(str_glue("atac: {atac_dsMin} - {atac_dsMax}."))
  query_data <- downsampleSeurat(seu = query_data, groupBy = query_cluster,
    minNum = query_downsample, maxNum = query_downsample)
}


# transfer annotation
RNA.anchors <- FindTransferAnchors(reference=ref_data, query=query_data, npcs=npc, dims=1:npc, features=NULL, k.anchor=kanchor, reduction=transfered, reference.reduction=ref_reduc)
predictions.Label <- TransferData(anchorset=RNA.anchors, refdata=ref_data@meta.data[,ref_cluster], query=query_data, dims=1:npc, weight.reduction=wet_reduc)
pt_query_umap_pred <- DimPlot(predictions.Label, reduction="umap", group.by="predicted.id", label=TRUE, label.size=4, repel=FALSE) + NoLegend()

# MapQuery
Ref.mapQuery <- MapQuery(anchorset=RNA.anchors, reference=ref_data, query=query_data,
    refdata=ref_data@meta.data[,ref_cluster], reference.reduction=ref_reduc, reduction.model="umap")
pt_ref_umap <- DimPlot(ref_data, reduction="umap", group.by=ref_cluster, label=TRUE, label.size=4, repel=FALSE) + NoLegend() + ggtitle("Reference annotations")
pt_qry_refumap <- DimPlot(Ref.mapQuery, reduction="ref.umap", group.by="predicted.id", label=TRUE, label.size=4, repel=FALSE) + NoLegend() + ggtitle("Query transferred labels")
pt_qry_umap <- DimPlot(Ref.mapQuery, reduction="umap", group.by="predicted.id", label=TRUE, label.size=4, repel=FALSE) + NoLegend() + ggtitle("Query transferred labels")

# overlap
#predictions.Label$celltype <- Idents(predictions.Label)
Label.table <- table(predictions.Label$predicted.id, predictions.Label@meta.data[,query_cluster])
Label.pct <- apply(Label.table, 2, function(x){x/sum(x)})


getSubclasstfScoreByclId <- function(scMat, top3cl, cl2subclass) {
 lapply(seq_along(top3cl), \(i) {
    l4 <- names(top3cl)[i]
    top3 <- top3cl[[i]]
    subclasses <- cl2subclass[names(top3), "subclass_id"] |>
      unique()
    scMat[as.character(subclasses), l4] |>
      setNames(object = _, nm = subclasses)
  }) |>
    setNames(object = _, nm = names(top3cl))
}

getPredictScoreMat <- function(tfquery) {
  predictAssay <- tfquery@assays$prediction.score.id
  score <- predictAssay@data
  return(score)
}

getPredictLabel <- function(tfquery) {
  labels <- tfquery$predicted.id
  scores <- tfquery$predicted.id.score
  cells <- colnames(tfquery)
  r <- data.frame(
    barcode = cells,
    predict = labels,
    score = scores
  )
  r$predict <- gsub("_", "-", r$predict) 
  rownames(r) <- cells
  return(r)
}

getMetaCol.Seurat <- function(seu, colnm) {
  ## use repeat colnm to get array
  ## instead of data.frame with one column
  g <- seu[[colnm]][[colnm]]
  names(g) <- colnames(seu)
  return(g)
}

toNamedArray.1dtable <- function(t) {
  tmp <- as.data.frame(t, stringsAsFactors = FALSE)
  r <- tmp$Freq
  names(r) <- tmp$Var1
  return(r)
}

Sys.setenv("_R_USE_PIPEBIND_" = TRUE)

getAvgVoteMat4SubclassByRefclId = function(seu, refCol = "ref_cluster", queryCol = "query_cluster", cl2subclass=NULL) {
      tfmat <- getPredictScoreMat(seu)
      
      refGroup <- rownames(tfmat)
#      scGroup <- unique(cl2subclass[refGroup, "subclass_id"])
      
      nrefCol <- nrow(tfmat)
      nquery <- ncol(tfmat)
      message(nrefCol, " groups in ref col: ", refCol)
      message(nquery, " cells are predicted.")
      
      tfLabels <- getPredictLabel(seu)
#      tfLabels$subclass_id <- cl2subclass[
#        tfLabels$predict, "subclass_id"]
      queryGroup <- getMetaCol.Seurat(seu, queryCol) |>
        x => x[tfLabels$barcode]
      ug <- unique(queryGroup)
      message("Find ", length(ug), " groups from ", queryCol)
      r <- future.apply::future_vapply(
        X = ug,
        FUN = \(g)  {
          r <- rep(0.0, length(refGroup))
          names(r) <- refGroup
          index <- queryGroup == g
          n <- sum(index)
          tmp <- table(tfLabels[index, "predict"]) |>
            toNamedArray.1dtable()
          r[names(tmp)] <- tmp / n
          return(r)
        },
        FUN.VALUE = rep(0.0, length(refGroup))
      )
      colnames(r) <- ug
      rownames(r) <- refGroup
      return(r)
    } ## end of fn

scMat <- getAvgVoteMat4SubclassByRefclId(predictions.Label, refCol=ref_cluster, queryCol=query_cluster)

#projdir <- here::here()
#rdir <- file.path(projdir, "scatac/scatac_mousebrain/package/R")
rdir <- file.path("/projects/ps-renlab2/kaw033/CEMBAv2/scatac/scatac_mousebrain/package/R")
import::from(.from = "ggtheme.R", .directory = rdir, dotplotTheme)

getMaxColScoreWithName <- function(mat) {
  maxScores <- apply(mat, 2, max)
  rnames <- rownames(mat)
  maxNames <- apply(mat, 2, \(x){
    rnames[which.max(x)]
  })
  names(maxScores) <- maxNames
  return(maxScores)
}
to3.matrix <- function(mat,
                       names = c("row", "column", "score"),
                       int2str = TRUE,
                       factor2str = TRUE) {
  r <- reshape2::melt(mat)
  colnames(r) <- names
  if (factor2str) {
    if (is.factor(r[,1])) {
      r[,1] <- as.character(r[,1])
    }
    if (is.factor(r[,2])) {
      r[,2] <- as.character(r[,2])
    }
  }
  if(int2str) {
    if(is.numeric(r[,1])) {
      r[,1] <- as.character(r[,1])
    }
    if(is.numeric(r[,2])) {
      r[,2] <- as.character(r[,2])
    }
  }
  return(r)
}

prepareDotPlot4TransferLabel <- function(tfmat,
                                         refOrder = NULL,
                                         names = c("row", "column", "score"),
                                         ignoreEmptyRef = TRUE) {
  maxScore <- getMaxColScoreWithName(mat = tfmat)
  query2ref <- data.frame(
    query = colnames(tfmat),
    ref = names(maxScore),
    row.names = colnames(tfmat)
  )
  if (ignoreEmptyRef) {
    message("remove refs not having query mapped to.")
    tfmat <- tfmat[rownames(tfmat) %in% query2ref$ref, ]
  }
  if (is.null(refOrder)) {
    message("refOrder is null, will use default numeric order for it.")
    refOrder <- rownames(tfmat)
    refOrder <- refOrder[order(as.integer(refOrder))]
  } else {
    refOrder <- refOrder[refOrder %in% rownames(tfmat)]
  }
  queryOrder <- query2ref$query[
    order(factor(query2ref$ref, levels = refOrder))]

  meltmat <- to3.matrix(tfmat, names)
  meltmat[,1] <- factor(meltmat[,1], levels = refOrder)
  meltmat[,2] <- factor(meltmat[,2], levels = queryOrder)
  # reduce size of meltmat
  meltmat <- meltmat[meltmat[,3] > 0, ]
  return(meltmat)
}

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
hclust_refs <- sort_hclust(hclust(dist(scMat), method = "ward.D2"))
#hclust_querys <- sort_hclust(hclust(dist(t(scMat)), method = "ward.D2"))

#mat.dotplot <- prepareDotPlot4TransferLabel(tfmat = scMat, ignoreEmptyRef = FALSE, names = c("ref", "query", "score"))
mat.dotplot <- prepareDotPlot4TransferLabel(tfmat = scMat, ignoreEmptyRef = FALSE, refOrder=hclust_refs$labels[hclust_refs$order], names = c("ref", "query", "score"))
mytheme <- dotplotTheme(legend.pos = "right")
lowSimScore<- quantile(mat.dotplot$score, 0.1)
#highSimScore <- max(mat.dotplot$score)
highSimScore <- quantile(mat.dotplot$score, 1)


palette <- "YlOrBr" #"Reds"
p.tfvote <- ggplot(data = mat.dotplot, aes(x = ref, y = query)) +
  geom_tile(aes(fill = score)) +
  mytheme +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size=8, color="black"),
    axis.text.y = element_text(size=8, color="black")) +
  #scale_fill_gradient(low = "white", high = "red",
  #  limits = c(lowSimScore, highSimScore), na.value = "white") +
  scale_fill_distiller(palette = palette, direction = 1, limits = c(lowSimScore, highSimScore), na.value = "white") + 
  #scale_size(range = c(0, 1)) +
  xlab(paste0(nrow(scMat), " clusters of ref")) +
  ylab(paste0(ncol(scMat), " clusters of query")) 

scMat <- scMat[match(levels(mat.dotplot$ref), rownames(scMat)), match(levels(mat.dotplot$query), colnames(scMat))]

find_unasigned <- function(scmat, matplot, cutoff=0.75) {
  unasigned_clusters <- names(which(colSums(scmat>=cutoff)==0))
  unasigned_pairs <- list()
  for (cluster in unasigned_clusters) {
    mat_cluster <- matplot[which(matplot$query == cluster),]
    match_clusters <- as.character(mat_cluster$ref[which(mat_cluster$score > 1/(nrow(mat_cluster)+1))])
    if (length(match_clusters)>1){
      unasigned_pairs[[cluster]] <- paste(cluster, match_clusters, sep=",")
    }
  }
  return(unasigned_pairs)
}

unasigned_pairs <- find_unasigned(scMat, mat.dotplot, cutoff=cutoff)

pdf(paste0(outF, ".TransferLabel.pdf"))
pt_query_umap_pred
pt_ref_umap
pt_qry_refumap
pt_qry_umap
pheatmap(Label.pct)
pheatmap(scMat, cluster_rows=FALSE, cluster_cols=FALSE)
p.tfvote
dev.off()

if (out_score) {
  write.csv(scMat, file=paste0(outF, ".TransferLabel.scMat.csv"), quote=FALSE)
  write.csv(mat.dotplot, file=paste0(outF, ".TransferLabel.matplot.csv"), quote=FALSE)
  unasigned_pairs_tbl <- as.data.frame(unlist(unasigned_pairs))
  write.csv(unasigned_pairs_tbl, file=paste0(outF, ".TransferLabel.unasigned.csv"), row.names=FALSE, col.names=FALSE, quote=FALSE)
}






