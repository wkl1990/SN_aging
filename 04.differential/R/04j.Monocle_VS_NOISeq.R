suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("monocle3"))

# compare monocle3 and NOISeq
Autochr.genetable_enzid <- read.csv(file="rds/RNA/preprocess/Autochr.genetable.entrezid.csv", row.names=1)
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)
RNA_subclass_forDEG <- intersect(rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell)>=100], rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0])

# scatter plot
subclass_label_id <- "327 Oligo NN"
print(subclass_label_id)
subclass_file <- gsub(" ", "_", subclass_label_id)
celltype_gene_fits <- readRDS(file=paste0("./rds/RNA/diff/after_integra/monocle3/subclass/", subclass_file, ".gene_fits.rds"))
celltype_fit_coefs <- coefficient_table(celltype_gene_fits)
celltype_age_terms <- celltype_fit_coefs %>% filter(term == "age")

print(paste0(subclass_label_id, " start!"))
# read deg
subclass_name <- gsub(" ", "_", subclass_label_id)
deg_file <- paste0("./rds/RNA/diff/after_integra/NOISeq/subclass/", subclass_name, "_noiseqbio_all.csv")
celltype_age_terms_noiseq <- read.csv(deg_file)


monocle3_NOISeq_join <- dplyr::inner_join(celltype_age_terms, celltype_age_terms_noiseq, by=join_by(name==X)) %>% mutate(theta=-theta, log2FC=-log2FC)

monocle_effect <- c("estimate",)
noiseq_effect <- c("log2FC")

plot_compare_monocle_noiseq <- function(join_df, monocle_name, noiseq_name, cor_method, title) {
  cor_val <- signif(cor(join_df[[monocle_name]], join_df[[noiseq_name]], method=cor_method), 2)
  pt <- ggplot(join_df, aes_string(x=monocle_name, y=noiseq_name)) +
    geom_point() + 
    geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + 
    annotate("text", x=min(join_df[[monocle_name]], na.rm=TRUE), y=max(join_df[[noiseq_name]], na.rm=TRUE), label=paste0(cor_method, " cor=", cor_val), size=5) + 
    xlab(paste0("Monocle3 differential expression statistics of age: ", monocle_name)) + 
    ylab(paste0("NOISeq differential expression statistics of age: ", noiseq_name)) + 
    ggtitle(title) + 
    theme_bw() + theme(axis.title=element_text(size=10, color="black"), axis.text=element_text(size=8, color="black"), plot.title=element_text(size=12, color="black", hjust=0.5))
  return(pt)
}

pt_monocle.estimate_noiseq.log2FC <- plot_compare_monocle_noiseq(monocle3_NOISeq_join, "estimate", "log2FC", "pearson", subclass_label_id)
ggsave(pt_monocle.estimate_noiseq.log2FC, file="./figures/fig2/DEG_monocle.estimate_noiseq.log2FC.scatter.pdf", width=8, height=8)

# count DEGs across subclasses 
# histogram plot
# monocle3-lm (q0.05)
monocle3_lm <- list()
monocle3_lm_q0.05_up <- list()
monocle3_lm_q0.05_down <- list()
for (subclass_label_id in RNA_subclass_forDEG) {
  subclass_name <- gsub(" ", "_", subclass_label_id)
  monocle3_lm[[subclass_label_id]] <- read.csv(paste0("/projects/ps-renlab2/kaw033/WK_aging/rds/RNA/diff/after_integra/monocle3/subclass/", subclass_name, ".age_terms_allgenes.csv"), row.names=1)
  monocle3_lm_q0.05_up[[subclass_label_id]] <- monocle3_lm[[subclass_label_id]]$name[intersect(which(monocle3_lm[[subclass_label_id]]$q_value<0.05), which(monocle3_lm[[subclass_label_id]]$estimate>0))]
  monocle3_lm_q0.05_down[[subclass_label_id]] <- monocle3_lm[[subclass_label_id]]$name[intersect(which(monocle3_lm[[subclass_label_id]]$q_value<0.05), which(monocle3_lm[[subclass_label_id]]$estimate<0))]
}
# NOISeq (prob0.95)
NOISeq <- list()
NOISeq_prob0.95_up <- list()
NOISeq_prob0.95_down <- list()
for (subclass_label_id in RNA_subclass_forDEG) {
  subclass_name <- gsub(" ", "_", subclass_label_id)
  NOISeq[[subclass_label_id]] <- read.csv(paste0("/projects/ps-renlab2/kaw033/WK_aging/rds/RNA/diff/after_integra/NOISeq/subclass/", subclass_name, "_noiseqbio_all.csv"))
  NOISeq_prob0.95_up[[subclass_label_id]] <- NOISeq[[subclass_label_id]]$X[intersect(which(NOISeq[[subclass_label_id]]$prob>0.95), which(NOISeq[[subclass_label_id]]$log2FC<0))]
  NOISeq_prob0.95_down[[subclass_label_id]] <- NOISeq[[subclass_label_id]]$X[intersect(which(NOISeq[[subclass_label_id]]$prob>0.95), which(NOISeq[[subclass_label_id]]$log2FC>0))]
}

# compare and plot
monocle3_lm_up <- do.call(c, monocle3_lm_q0.05_up) %>% as.data.frame %>% rename("."="name") %>% mutate(gene=paste0(name, "_up")) 
monocle3_lm_down <- do.call(c, monocle3_lm_q0.05_down) %>% as.data.frame %>% rename("."="name") %>% mutate(gene=paste0(name, "_up")) 
monocle3_lm_combined <- rbind(monocle3_lm_up, monocle3_lm_down)

monocle3_lm_combined %>% group_by(gene) %>% summarise(count=n()) %>% arrange(desc(count)) -> monocle3_lm_combined_count
monocle3_lm_combined %>% group_by(name) %>% summarise(count=n()) %>% arrange(desc(count)) -> monocle3_lm_combined_count1

monocle3_lm_topgene <- monocle3_lm_combined_count %>% filter(count==max(count)) %>% pull(gene)

pt_monocle_DEGcount <- ggplot(monocle3_lm_combined_count, aes(x=count)) +
  geom_histogram() + scale_x_continuous(breaks=seq(1, max(monocle3_lm_combined_count$count), by=1)) + 
  annotate("text", x=max(monocle3_lm_combined_count$count), y=nrow(monocle3_lm_combined_count)/10, label=paste(monocle3_lm_topgene, collapse="\n"), size=4, color="black") +
  xlab("Number of cell subclasses") + ylab("Count") + theme_classic() +
  theme(axis.title=element_text(size=10, color="black"), axis.text=element_text(size=8, color="black"), plot.title=element_text(size=12, color="black", hjust=0.5)) 
ggsave(pt_monocle_DEGcount, file="./figures/fig2/DEG_monocle_DEGcount.histogram.pdf", width=8, height=8)


NOISeq_up <- do.call(c, NOISeq_prob0.95_up) %>% as.data.frame %>% rename("."="name") %>% mutate(gene=paste0(name, "_up")) 
NOISeq_down <- do.call(c, NOISeq_prob0.95_down) %>% as.data.frame %>% rename("."="name") %>% mutate(gene=paste0(name, "_down")) 
NOISeq_combined <- rbind(NOISeq_up, NOISeq_down)

NOISeq_combined %>% group_by(gene) %>% summarise(count=n()) %>% arrange(desc(count)) -> NOISeq_combined_count
NOISeq_combined %>% group_by(name) %>% summarise(count=n()) %>% arrange(desc(count)) -> NOISeq_combined_count1

NOISeq_topgene <- NOISeq_combined_count %>% filter(count==max(count)) %>% pull(gene)

pt_NOISeq_DEGcount <- ggplot(NOISeq_combined_count, aes(x=count)) +
  geom_histogram() + scale_x_continuous(breaks=seq(1, max(NOISeq_combined_count$count), by=1)) +
  annotate("text", x=max(NOISeq_combined_count$count), y=nrow(NOISeq_combined_count)/10, label=paste(NOISeq_topgene, collapse="\n"), size=4, color="black") +
  xlab("Number of cell subclasses") + ylab("Count") + theme_classic() + 
  theme(axis.title=element_text(size=10, color="black"), axis.text=element_text(size=8, color="black"), plot.title=element_text(size=12, color="black", hjust=0.5))
ggsave(pt_NOISeq_DEGcount, file="./figures/fig2/DEG_NOISeq_DEGcount.histogram.pdf", width=8, height=8)
