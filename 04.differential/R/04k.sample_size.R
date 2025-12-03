suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("monocle3"))
suppressPackageStartupMessages(library("scPS"))

# sample size
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)
RNA$age <- case_match(RNA$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18)
RNA$rep <- case_match(RNA$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2")

# using all data
counts <- RNA@assays$RNA$counts
cell.info <- RNA@meta.data %>% mutate(age=factor(age, levels=c(2,6,12,18)), subclass=subclass_label_id) %>% select(sampleID, age, subclass) 


# test
RNA_subclass_sample_cell <- table(RNA$subclass_label_id, RNA$sampleID)
RNA_subclass_forDEG <- intersect(rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell)>=100], rownames(RNA_subclass_sample_cell)[rowSums(RNA_subclass_sample_cell<10)==0])

geneObject <- estPreParas.multi(counts, cell.info,
                                id="sampleID", x1="age", cellcluster="subclass",
                                cells.interesting="327 Oligo NN")

Genes.tested <- geneCandidate(geneObject, thres.nonZeroPs = 0.1, nGenesCandidate = 1000, propDEGs = 0.01)
mean1 <- Genes.tested$`327 Oligo NN`$mean.control
icc <- Genes.tested$`327 Oligo NN`$icc
hf <- Genes.tested$`327 Oligo NN`$hf


dir.create("./figures/fig2/samplesize/test_params/", showWarnings = FALSE)
for (sub in c(2,4,8)) {
	print(paste0("Sample size calculation for sub=", sub))
	mean1 <- Genes.tested$`327 Oligo NN`$mean.control
	icc <- Genes.tested$`327 Oligo NN`$icc
	hf <- Genes.tested$`327 Oligo NN`$hf
	esizes <- seq(2.0, 2.5, 0.1)
	list3 <- lapply(esizes, function(x) {
	  FC <- c(rep(x, 50), rep(1, 950))
	  size.view <- sizeCal(low.up.m=c(sub,sub), low.up.n=c(100,1000), ePower=0.8, FDR=0.05,
	                        grid.m=1, grid.n=100, r=1, rc=1, total=NULL,
	                        vvmean1=mean1, FC=FC, vvrho=icc, hf=hf)
	  cbind(x=x, size.view$m.n.power)
	})
	dat2 <- do.call(rbind, list3); ePower <- 0.8
	tryCatch({
	  fig <- ggplot(dat2, aes(x=x, y=n, fill=power)) +
		  geom_point(size=10, shape=21, colour = "transparent") +
		  geom_text(aes(label = round(power, 2), color = ifelse(power > ePower, "blue", "red"), fontface=2),
		            size = 3.2, show.legend = FALSE) +
		  scale_color_manual(values = c("blue", "red")) +
		  scale_fill_gradient(low = "yellow", high = "green") +
		  scale_x_continuous(breaks = dat2$x) +
		  scale_y_continuous(breaks = dat2$n) +
		  xlab("Effect size (FC)") +
		  ylab("No. of cells per stage") +
		  theme_minimal()
	  ggsave(fig, file=paste0("./figures/fig2/samplesize/test_params/sample_size_power_sub", sub, ".pdf"), width=8, height=8)
	}, error=function(e){print(paste0("Error in plotting for sub=", sub))})

}
