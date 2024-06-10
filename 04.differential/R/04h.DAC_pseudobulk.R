suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))

# read data
input <- "rds/ATAC/after_integra/seurat/SN_integra.final.peak.srt.pmat.seurat4.rds"
ATAC <- readRDS(input)
Idents(ATAC) <- ATAC@meta.data$subclass_label_id
ATAC$sampleID <- factor(ATAC$orig.ident, levels=c("2m_rep1", "2m_rep2", "6m_rep1", "6m_rep2", "12m_rep1", "12m_rep2", "18m_rep1", "18m_rep2"))

ATAC$age <- factor(case_match(ATAC$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18), levels=c(2, 6, 12, 18))
ATAC$rep <- factor(case_match(ATAC$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2"), levels=c("rep1", "rep2"))
ATAC$age_rep <- factor(case_match(ATAC$sampleID, "2m_rep1" ~ 2.1, "6m_rep1" ~ 6.1, "12m_rep1" ~ 12.1, "18m_rep1" ~ 18.1, "2m_rep2" ~ 2.2, "6m_rep2" ~ 6.2, "12m_rep2" ~ 12.2, "18m_rep2" ~ 18.2), levels=c(2.1,2.2,6.1,6.2,12.1,12.2,18.1,18.2))

# seurat pseudo-bulk
pseudo_bulk <- function(data=ATAC, celltype=celltype, idents="celltype", assay="ATAC", slot="count", group.by="sampleID"){
  Idents(data) <- ATAC@meta.data[[idents]]
  ATAC_celltype <- subset(data, idents=celltype)
  ATAC_celltype_pseudobulk <- AggregateExpression(ATAC_celltype, assays=assay, slot=slot, group.by=group.by)$ATAC
  colnames(ATAC_celltype_pseudobulk) <- gsub("^g", "", colnames(ATAC_celltype_pseudobulk))
  colnames(ATAC_celltype_pseudobulk) <- gsub("-", "_", colnames(ATAC_celltype_pseudobulk))
  return(ATAC_celltype_pseudobulk)
}

pseudo_filter <- function(data=ATAC_pseudobulk, depth=1000000, count=1, samples=1){
  data <- data[rowSums(data >= count) >= samples, ]
  data <- data[, colSums(data) >= depth, drop=FALSE]
  return(data)
}


# meta table
meta_table <- function(data=ATAC, celltype=celltype, sample="sampleID", ident="celltype", select=c("sampleID", "age", "rep", "age_rep", "sample_cells")){
  Idents(data) <- ATAC@meta.data[[ident]]
  ATAC_celltype <- subset(data, idents=celltype)
  sample_cells <- table(data@meta.data[,sample]) %>%  as.vector()
  names(sample_cells) <- names(table(data@meta.data[,sample]))
  m <- match(names(sample_cells), data@meta.data[,sample])
  ATAC_celltype_meta <- data.frame(ATAC@meta.data[m, ], sample_cells, row.names = NULL) %>% dplyr::select(select)
  t <- table(ATAC_celltype@meta.data[,sample], ATAC_celltype@meta.data[,ident])
  cell_counts <- t[, which(colnames(t) == celltype)]
  ATAC_celltype_meta$n_cells <- cell_counts[match(ATAC_celltype_meta$sampleID, names(cell_counts))]
  sumdata <- ATAC_celltype@meta.data %>% group_by(sampleID) %>% summarise(nCount=sum(nCount_ATAC), nFeature=sum(nFeature_ATAC))
  ATAC_celltype_meta <- plyr::join(ATAC_celltype_meta, sumdata, by = intersect(names(ATAC_celltype_meta), names(sumdata)))
  rownames(ATAC_celltype_meta) <- ATAC_celltype_meta[,sample]
  return(ATAC_celltype_meta)
}


ATAC_subclass_sample_cell <- table(ATAC$subclass_label_id, ATAC$sampleID)
ATAC_subclass_forDEG <- intersect(rownames(ATAC_subclass_sample_cell)[rowSums(ATAC_subclass_sample_cell)>=100], rownames(ATAC_subclass_sample_cell)[rowSums(ATAC_subclass_sample_cell<10)==0])

# Run the script on all clusters comparing 18m relative to 2m (filter sex chromosomes)
contrast_samples <- c("2m_rep1", "2m_rep2", "18m_rep1", "18m_rep2")

for (subclass_label_id in ATAC_subclass_forDEG){
  print(subclass_label_id)
  name <- gsub(" ", "_", subclass_label_id)
  ATAC_celltype_pseudobulk <- pseudo_bulk(data=ATAC, celltype=subclass_label_id, ident="subclass_label_id", assay="ATAC", slot="count", group.by="sampleID")
  ATAC_celltype_pseudobulk_filter <- pseudo_filter(data=ATAC_celltype_pseudobulk, depth=1000, count=1, samples=1)
  write.csv(ATAC_celltype_pseudobulk_filter, file=paste0("rds/ATAC/after_integra/diff/pseudobulk/", name, "_pseudobulk.csv"), quote = FALSE)
  ATAC_celltype_pseudobulk_filter_autochr <- ATAC_celltype_pseudobulk_filter[grep("chr[1-9]", rownames(ATAC_celltype_pseudobulk_filter)), ]
  write.csv(ATAC_celltype_pseudobulk_filter_autochr, file=paste0("rds/ATAC/after_integra/diff/pseudobulk/", name, "_pseudobulk_autochr.csv"), quote = FALSE)
  ATAC_celltype_meta <- meta_table(data=ATAC, celltype=subclass_label_id, sample="sampleID", ident="subclass_label_id", select=c("sampleID", "age", "rep", "age_rep", "sample_cells"))
  write.csv(ATAC_celltype_meta, file=paste0("rds/ATAC/after_integra/diff/pseudobulk/", name, "_meta.csv"), quote = FALSE)
  ATAC_celltype_pseudobulk_filter_contrast <- ATAC_celltype_pseudobulk_filter[, colnames(ATAC_celltype_pseudobulk_filter) %in% contrast_samples]
  write.csv(ATAC_celltype_pseudobulk_filter_contrast, file=paste0("rds/ATAC/after_integra/diff/pseudobulk/", name, "_pseudobulk_contrast.csv"), quote = FALSE)
  ATAC_celltype_pseudobulk_filter_contrast_autochr <- ATAC_celltype_pseudobulk_filter_contrast[grep("chr[1-9]", rownames(ATAC_celltype_pseudobulk_filter_contrast)), ]
  write.csv(ATAC_celltype_pseudobulk_filter_contrast_autochr, file=paste0("rds/ATAC/after_integra/diff/pseudobulk/", name, "_pseudobulk_contrast_autochr.csv"), quote = FALSE)
  ATAC_celltype_meta_contrast <- ATAC_celltype_meta[ATAC_celltype_meta$sampleID %in% contrast_samples, ]  
  write.csv(ATAC_celltype_meta_contrast, file=paste0("rds/ATAC/after_integra/diff/pseudobulk/", name, "_meta_contrast.csv"), quote = FALSE)
  print(paste0(subclass_label_id, " finished!"))
}


# try NOISeq
library(NOISeq)
library(EDASeq)

ATAC_c8_pseudobulk <- pseudo_bulk(data=ATAC, celltype=subclass_label_id, ident="subclass_label_id", assay="ATAC", slot="count", group.by="sampleID")
ATAC_c8_pseudobulk_filter <- pseudo_filter(data=ATAC_c8_pseudobulk, depth=1000, count=1, samples=1)

genes_annot <- data.frame(gene=rownames(ATAC_c8_pseudobulk))
rownames(genes_annot) <- genes_annot[,1]
genes_annot %>% tidyr::separate(gene, sep="[:-]", into=c("chr", "start", "end"), remove=FALSE) -> genes_annot

contrast_samples <- c("2m_rep1", "2m_rep2", "18m_rep1", "18m_rep2")

noiseqbio_deg <- function(data, meta, annot, gene="gene", chr="chr", start="start", end="end", factor="age", name="celltype", path){
  celltypedata <- NOISeq::readData(data=data, chromosome=genes_annot[,c(chr,start,end)], factors=meta)
  QCreport(celltypedata, file=paste0(path, name, "_QCreport.pdf"), samples = NULL, factor = factor, norm = FALSE)
  celltypenoiseqbio = noiseqbio(celltypedata, k = 0.5, norm = "rpkm", factor = factor, lc = 1, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345, filter = 1)
  celltypenoiseqbio.deg = degenes(celltypenoiseqbio, q = 0.95, M = NULL)
  celltypenoiseqbio.all = degenes(celltypenoiseqbio, q = 0, M = NULL)
  celltypenoiseqbio.degDown = degenes(celltypenoiseqbio, q = 0.95, M = "up")
  celltypenoiseqbio.degUp = degenes(celltypenoiseqbio, q = 0.95, M = "down")
  write.csv(celltypenoiseqbio.all, file=paste0(path, name, "_noiseqbio_all.csv"), quote = FALSE)
  write.csv(celltypenoiseqbio.deg, file=paste0(path, name, "_noiseqbio_deg.csv"), quote = FALSE)
  write.csv(celltypenoiseqbio.degDown, file=paste0(path, name, "_noiseqbio_degDown.csv"), quote = FALSE)
  write.csv(celltypenoiseqbio.degUp, file=paste0(path, name, "_noiseqbio_degUp.csv"), quote = FALSE)
  pdf(paste0(path, name, "_DEGplot.pdf"))
  DE.plot(celltypenoiseqbio, q = 0.95, graphic = "expr", log.scale = TRUE)
  DE.plot(celltypenoiseqbio, q = 0.95, graphic = "MD")
  #DE.plot(celltypenoiseqbio, chromosomes = c(1, 2), log.scale = TRUE, join = FALSE, q = 0.95, graphic = "chrom")
  DE.plot(celltypenoiseqbio, chromosomes = NULL, q = 0.95, graphic = "distr")
  dev.off()
}

# DEG (filter sex chromosome)
for (subclass_label_id in ATAC_subclass_forDEG){
  print(subclass_label_id)
  name <- gsub(" ", "_", subclass_label_id)
  ATAC_celltype_pseudobulk_filter_contrast <- read.csv(paste0("rds/ATAC/after_integra/diff/pseudobulk/", name, "_pseudobulk_contrast.csv"), row.names=1)
  ATAC_celltype_meta_contrast <- read.csv(paste0("rds/ATAC/after_integra/diff/pseudobulk/", name, "_meta_contrast.csv"), row.names=1)
  #ATAC_celltype_pseudobulk_filter_contrast_filteATAC <- ATAC_celltype_pseudobulk_filter_contrast[!rownames(ATAC_celltype_pseudobulk_filter_contrast) %in% lengthNA_genes, ]
  ATAC_celltype_pseudobulk_filter_contrast_autochr <- ATAC_celltype_pseudobulk_filter_contrast[grep("chr[1-9]", rownames(ATAC_celltype_pseudobulk_filter_contrast)), ]
  ATAC_celltype_pseudobulk_filter_contrast_autochr_filter <- pseudo_filter(data=ATAC_celltype_pseudobulk_filter_contrast_autochr, depth=1000, count=1, samples=1)
  message("!!!!!", subclass_label_id, " expressed peaks number: ", nrow(ATAC_celltype_pseudobulk_filter_contrast_autochr), "!!!!!")
  message("!!!!!", subclass_label_id, " expressed peaks number: ", nrow(ATAC_celltype_pseudobulk_filter_contrast_autochr_filter), "!!!!!")
  # filter subclass peaks
  celltype_peaks <- read.table(paste0("/peak_calling/after_integra/final/subclass.final.peak.srt/", name, ".bed"))
  celltype_called_peaks <- paste0(celltype_peaks$V1,":",celltype_peaks$V2,"-",celltype_peaks$V3)
  message("!!!!!", subclass_label_id, " called peaks number: ", length(celltype_called_peaks), "!!!!!")
  ATAC_celltype_pseudobulk_filter_contrast_autochr <- ATAC_celltype_pseudobulk_filter_contrast_autochr[which(rownames(ATAC_celltype_pseudobulk_filter_contrast_autochr) %in% celltype_called_peaks), ]
  message("!!!!!", subclass_label_id, " expressed peaks number: ", nrow(ATAC_celltype_pseudobulk_filter_contrast_autochr), "!!!!!")  
  tryCatch({noiseqbio_deg(data=ATAC_celltype_pseudobulk_filter_contrast_autochr, meta=ATAC_celltype_meta_contrast, annot=genes_annot, 
    gene="gene", chr="chr", start="start", end="end", factor="age", name=name, path="rds/ATAC/after_integra/diff/NOISeq/")}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  print(paste0(subclass_label_id, " finished!"))
}


NOISeq_DEGs_statistics <- function(file, var="prob", cutoff=1, direction="up"){
  DEGs_tbl <- read.csv(file)
  if (direction == "up"){
    n_sig_up <- sum(DEGs_tbl[,var] > cutoff & DEGs_tbl[,"log2FC"] < 0, na.rm=TRUE)
    return(n_sig_up)
  } else if (direction == "down"){
    n_sig_down <- sum(DEGs_tbl[,var] > cutoff & DEGs_tbl[,"log2FC"] > 0, na.rm=TRUE)
    return(n_sig_down)
  } else {
    n_sig <- sum(DEGs_tbl[,var] > cutoff, na.rm=TRUE)
    return(n_sig)
  }
}

DEGs_NOISeq_prob0.95_up <- c()
DEGs_NOISeq_prob0.95_down <- c()
for (subclass_label_id in ATAC_subclass_forDEG) {
  name <- gsub(" ", "_", subclass_label_id)
  file <- paste0("rds/ATAC/after_integra/diff/NOISeq/", name, "_noiseqbio_all.csv")
  n_prob0.95_up <- NOISeq_DEGs_statistics(file, var="prob", cutoff=0.95, direction="up")
  n_prob0.95_down <- NOISeq_DEGs_statistics(file, var="prob", cutoff=0.95, direction="down")
  DEGs_NOISeq_prob0.95_up <- c(DEGs_NOISeq_prob0.95_up, n_prob0.95_up)
  DEGs_NOISeq_prob0.95_down <- c(DEGs_NOISeq_prob0.95_down, n_prob0.95_down)
}
DEGs_NOISeq <- cbind(DEGs_NOISeq_prob0.95_up, DEGs_NOISeq_prob0.95_down)
rownames(DEGs_NOISeq) <- ATAC_subclass_forDEG
colnames(DEGs_NOISeq) <- c("prob0.95_up", "prob0.95_down")
DEGs_NOISeq_long <- melt(DEGs_NOISeq) %>% data.frame() %>% dplyr::rename("celltype"="Var1", "class"="Var2") %>% 
  dplyr::mutate(class=factor(class, levels=c("prob0.95_up", "prob0.95_down"))) %>% 
  dplyr::mutate(celltype=factor(celltype, levels=ATAC_subclass_forDEG))

pt_NOISeqDEG <- ggplot(DEGs_NOISeq_long, aes(fill=class, y=value, x=celltype)) + geom_bar(position="dodge", stat="identity", width=0.5) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) 
pt_NOISeqDEG_theme <- pt_NOISeqDEG + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="number of DEGs", fill="cutoff") + 
  ggtitle("Number of DEGs per cell type") + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette="Paired")
ggsave(pt_NOISeqDEG_theme, file="rds/ATAC/after_integra/diff/NOISeq/pt_NOISeqDEG.pdf", width=16, height=4)

