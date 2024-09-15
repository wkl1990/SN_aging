suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))

# read data
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
RNA$age <- factor(case_match(RNA$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18), levels=c(2, 6, 12, 18))
RNA$rep <- factor(case_match(RNA$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2"), levels=c("rep1", "rep2"))
RNA$age_rep <- factor(case_match(RNA$sampleID, "2m_rep1" ~ 2.1, "6m_rep1" ~ 6.1, "12m_rep1" ~ 12.1, "18m_rep1" ~ 18.1, "2m_rep2" ~ 2.2, "6m_rep2" ~ 6.2, "12m_rep2" ~ 12.2, "18m_rep2" ~ 18.2), levels=c(2.1,2.2,6.1,6.2,12.1,12.2,18.1,18.2))

# seurat pseudo-bulk
pseudo_bulk <- function(data=RNA, celltype=celltype, idents="celltype", assay="RNA", slot="count", group.by="sampleID"){
  Idents(data) <- data@meta.data[[idents]]
  RNA_celltype <- subset(data, idents=celltype)
  RNA_celltype_pseudobulk <- AggregateExpression(RNA_celltype, assays=assay, slot=slot, group.by=group.by)$RNA
  colnames(RNA_celltype_pseudobulk) <- gsub("^g", "", colnames(RNA_celltype_pseudobulk))
  colnames(RNA_celltype_pseudobulk) <- gsub("-", "_", colnames(RNA_celltype_pseudobulk))
  return(RNA_celltype_pseudobulk)
}

pseudo_filter <- function(data=RNA_pseudobulk, depth=1000000, count=1, samples=1){
  data <- data[rowSums(data >= count) >= samples, ]
  data <- data[, colSums(data) >= depth, drop=FALSE]
  return(data)
}


# meta table
meta_table <- function(data=RNA, celltype=celltype, sample="sampleID", ident="celltype", select=c("sampleID", "age", "rep", "age_rep", "sample_cells")){
  Idents(data) <- data@meta.data[[ident]]
  RNA_celltype <- subset(data, idents=celltype)
  sample_cells <- table(data@meta.data[,sample]) %>%  as.vector()
  names(sample_cells) <- names(table(data@meta.data[,sample]))
  m <- match(names(sample_cells), data@meta.data[,sample])
  RNA_celltype_meta <- data.frame(data@meta.data[m, ], sample_cells, row.names = NULL) %>% dplyr::select(select)
  t <- table(RNA_celltype@meta.data[,sample], RNA_celltype@meta.data[,ident])
  cell_counts <- t[, which(colnames(t) == celltype)]
  RNA_celltype_meta$n_cells <- cell_counts[match(RNA_celltype_meta$sampleID, names(cell_counts))]
  sumdata <- RNA_celltype@meta.data %>% group_by(sampleID) %>% summarise(nCount=sum(nCount_RNA), nFeature=sum(nFeature_RNA), mt_perc=mean(percent.mt))
  RNA_celltype_meta <- plyr::join(RNA_celltype_meta, sumdata, by = intersect(names(RNA_celltype_meta), names(sumdata)))
  rownames(RNA_celltype_meta) <- RNA_celltype_meta[,sample]
  return(RNA_celltype_meta)
}


# Run the script on all clusters comparing 18m relative to 2m (filter sex chromosomes)
Autochr.genetable <- read.csv(file="rds/RNA/preprocess/Autochr.genetable.csv")
contrast_samples <- c("2m_rep1", "2m_rep2", "18m_rep1", "18m_rep2")


#RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)[,contrast_samples]
#RNA_subclass_forDEG <- rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0]
RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)
RNA_subclass_forDEG <- intersect(rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell)>=100], rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0])


for (subclass_label_id in RNA_subclass_forDEG){
  print(subclass_label_id)
  name <- gsub(" ", "_", subclass_label_id)
  RNA_subclass_pseudobulk <- pseudo_bulk(data=RNA, celltype=subclass_label_id, idents="subclass_label_id", assay="RNA", slot="count", group.by="sampleID")
  RNA_subclass_pseudobulk_filter <- pseudo_filter(data=RNA_subclass_pseudobulk, depth=1000, count=1, samples=1)
  write.csv(RNA_subclass_pseudobulk_filter, file=paste0("rds/RNA/diff/after_integra/pseudo_bulk/subclass/", name, "_pseudobulk.csv"), quote = FALSE)
  RNA_subclass_pseudobulk_filter_autochr <- RNA_subclass_pseudobulk_filter[rownames(RNA_subclass_pseudobulk_filter) %in% Autochr.genetable$gene_name, ]
  write.csv(RNA_subclass_pseudobulk_filter_autochr, file=paste0("rds/RNA/diff/after_integra/pseudo_bulk/subclass/", name, "_pseudobulk_autochr.csv"), quote = FALSE)
  RNA_subclass_meta <- meta_table(data=RNA, celltype=subclass_label_id, sample="sampleID", ident="subclass_label_id", select=c("sampleID", "age", "rep", "age_rep", "sample_cells"))
  write.csv(RNA_subclass_meta, file=paste0("rds/RNA/diff/after_integra/pseudo_bulk/subclass/", name, "_meta.csv"), quote = FALSE)
  RNA_subclass_pseudobulk_filter_contrast <- RNA_subclass_pseudobulk_filter[, colnames(RNA_subclass_pseudobulk_filter) %in% contrast_samples]
  write.csv(RNA_subclass_pseudobulk_filter_contrast, file=paste0("rds/RNA/diff/after_integra/pseudo_bulk/subclass/", name, "_pseudobulk_contrast.csv"), quote = FALSE)
  RNA_subclass_pseudobulk_filter_contrast_autochr <- RNA_subclass_pseudobulk_filter_contrast[rownames(RNA_subclass_pseudobulk_filter_contrast) %in% Autochr.genetable$gene_name, ]
  write.csv(RNA_subclass_pseudobulk_filter_contrast_autochr, file=paste0("rds/RNA/diff/after_integra/pseudo_bulk/subclass/", name, "_pseudobulk_contrast_autochr.csv"), quote = FALSE)
  RNA_subclass_meta_contrast <- RNA_subclass_meta[RNA_subclass_meta$sampleID %in% contrast_samples, ] 
  write.csv(RNA_subclass_meta_contrast, file=paste0("rds/RNA/diff/after_integra/pseudo_bulk/subclass/", name, "_meta_contrast.csv"), quote = FALSE)
#  tryCatch({get_dds_resultsAvsB(data=RNA_celltype_pseudobulk_filter_contrast_autochr, meta=RNA_celltype_meta_contrast, A="18", B="2", padj_cutoff=0.05, name=celltype)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#  tryCatch({get_dds_LRTresults(data=RNA_celltype_pseudobulk_filter_autochr, meta=RNA_celltype_meta, name=celltype)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  print(paste0(subclass_label_id, " finished!"))
}

# try NOISeq
library(NOISeq)
library(EDASeq)

genes_annot <- read.table("annot/genes.txt")
rownames(genes_annot) <- genes_annot[,1]
colnames(genes_annot) <- c("gene", "chr", "start", "end", "strand", "gene_type", "gene_name")
genes_annot$dup[duplicated(genes_annot$gene_name)] <- 1
genes_annot <- genes_annot %>% mutate(gene_name_mod=case_match(dup, 1~paste0(gene_name, ".", dup), .default=gene_name))
genes_annot$dup[duplicated(genes_annot$gene_name_mod)] <- 2
genes_annot <- genes_annot %>% mutate(gene_name_mod=case_match(dup, 2~paste0(gene_name, ".", dup), .default=gene_name_mod))

#genes_annot_len <- getGeneLengthAndGCContent(genes_annot$gene, org="mmu", mode="org.db")
#genes_annot_len_biomart <- getGeneLengthAndGCContent(genes_annot$gene, org="mmu")
#genes_annot_filter <- genes_annot %>% filter(gene_type %in% c("protein_coding", "lincRNA"))
genes_annot_len <- read.table("annot/genes.gtf.gene_length.txt", header=TRUE)
genes_annot_len_biomart <- getGeneLengthAndGCContent(genes_annot$gene, org="mmu")
genes_annot$gc <- genes_annot_len_biomart[match(genes_annot$gene, rownames(genes_annot_len_biomart)), "gc"]
genes_annot$length_gtftools <- genes_annot_len$longest_isoform[match(genes_annot$gene, genes_annot_len$gene)]
genes_annot$length_biomart <- genes_annot_len_biomart[match(genes_annot$gene, rownames(genes_annot_len_biomart)), "length"]
genes_annot <- genes_annot %>% mutate(length=case_when(is.na(length_gtftools)~length_biomart, TRUE~length_gtftools))
rownames(genes_annot) <- genes_annot$gene_name_mod
lengthNA_genes <- genes_annot$gene_name_mod[is.na(genes_annot$length)]

contrast_samples <- c("2m_rep1", "2m_rep2", "18m_rep1", "18m_rep2")

noiseqbio_deg <- function(data, meta, annot, gene="gene_name_mod", len="length", gc="gc", biotype="gene_type", chr="chr", start="start", end="end", factor="age", name="celltype", path){
  celltypedata <- NOISeq::readData(data=data, length=annot[,c(gene,len)], gc=genes_annot[,c(gene,gc)], biotype=genes_annot[,c(gene,biotype)], chromosome=genes_annot[,c(chr,start,end)], factors=meta)
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
Autochr.genetable <- read.csv(file="rds/RNA/preprocess/Autochr.genetable.csv")

for (subclass_label_id in RNA_subclass_forDEG){
  print(subclass_label_id)
  name <- gsub(" ", "_", subclass_label_id)
  RNA_subclass_pseudobulk_filter_contrast <- read.csv(paste0("rds/RNA/diff/after_integra/pseudo_bulk/subclass/", name, "_pseudobulk_contrast.csv"), row.names=1)
  RNA_subclass_meta_contrast <- read.csv(paste0("rds/RNA/diff/after_integra/pseudo_bulk/subclass/", name, "_meta_contrast.csv"), row.names=1)
  RNA_subclass_pseudobulk_filter_contrast_filterNA <- RNA_subclass_pseudobulk_filter_contrast[!rownames(RNA_subclass_pseudobulk_filter_contrast) %in% lengthNA_genes, ]
  RNA_subclass_pseudobulk_filter_contrast_filterNA_autochr <- RNA_subclass_pseudobulk_filter_contrast_filterNA[rownames(RNA_subclass_pseudobulk_filter_contrast_filterNA) %in% Autochr.genetable$gene_name, ]
  tryCatch({noiseqbio_deg(data=RNA_subclass_pseudobulk_filter_contrast_filterNA_autochr, meta=RNA_subclass_meta_contrast, annot=genes_annot, 
    gene="gene_name_mod", len="length", gc="gc", biotype="gene_type", chr="chr", start="start", end="end", factor="age", name=name, path="rds/RNA/diff/after_integra/NOISeq/subclass/")}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
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

subclass_DEGs_NOISeq_prob0.95_up <- c()
subclass_DEGs_NOISeq_prob0.95_down <- c()
for (subclass_label_id in RNA_subclass_forDEG) {
  name <- gsub(" ", "_", subclass_label_id)
  file <- paste0("rds/RNA/diff/after_integra/NOISeq/subclass/", name, "_noiseqbio_all.csv")
  n_prob0.95_up <- NOISeq_DEGs_statistics(file, var="prob", cutoff=0.95, direction="up")
  n_prob0.95_down <- NOISeq_DEGs_statistics(file, var="prob", cutoff=0.95, direction="down")
  subclass_DEGs_NOISeq_prob0.95_up <- c(subclass_DEGs_NOISeq_prob0.95_up, n_prob0.95_up)
  subclass_DEGs_NOISeq_prob0.95_down <- c(subclass_DEGs_NOISeq_prob0.95_down, n_prob0.95_down)
}
subclass_DEGs_NOISeq <- cbind(subclass_DEGs_NOISeq_prob0.95_up, subclass_DEGs_NOISeq_prob0.95_down)
rownames(subclass_DEGs_NOISeq) <- RNA_subclass_forDEG
colnames(subclass_DEGs_NOISeq) <- c("prob0.95_up", "prob0.95_down")
subclass_DEGs_NOISeq_long <- melt(subclass_DEGs_NOISeq) %>% data.frame() %>% dplyr::rename("cellsubclass"="Var1", "class"="Var2") %>% 
  dplyr::mutate(class=factor(class, levels=c("prob0.95_up", "prob0.95_down"))) %>% 
  dplyr::mutate(cellsubclass=factor(cellsubclass, levels=RNA_subclass_forDEG))

pt_subclass_NOISeqDEG <- ggplot(subclass_DEGs_NOISeq_long, aes(fill=class, y=value, x=cellsubclass)) + geom_bar(position="dodge", stat="identity", width=0.5) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) 
pt_subclass_NOISeqDEG_theme <- pt_subclass_NOISeqDEG + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="number of DEGs", fill="cutoff") + 
  ggtitle("Number of DEGs per cell subclass") + theme(plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette="Paired")


