#!/usr/bin/env Rscript

# read data
input <- "./rds/RNA/cluster/redoround2/RNA.combined.allen.integration.final.rds"
RNA <- readRDS(input)

RNA$age <- case_match(RNA$sampleID, c("2m_rep1", "2m_rep2") ~ 2, c("6m_rep1", "6m_rep2") ~ 6, c("12m_rep1", "12m_rep2") ~ 12, c("18m_rep1", "18m_rep2") ~ 18)
RNA$rep <- case_match(RNA$sampleID, c("2m_rep1", "6m_rep1", "12m_rep1", "18m_rep1") ~ "rep1", c("2m_rep2", "6m_rep2", "12m_rep2", "18m_rep2") ~ "rep2")
RNA$age_rep <- case_match(RNA$sampleID, "2m_rep1" ~ 2.1, "6m_rep1" ~ 6.1, "12m_rep1" ~ 12.1, "18m_rep1" ~ 18.1, "2m_rep2" ~ 2.2, "6m_rep2" ~ 6.2, "12m_rep2" ~ 12.2, "18m_rep2" ~ 18.2)

# cell proportion
RNA@meta.data %>% select(sampleID, class_label_id, subclass_label_id, rep, age, age_rep) -> RNA.celldict
RNA@meta.data %>% select(subclass_label_id, subclass_color) %>% distinct %>% arrange(subclass_label_id) -> RNA_subclass_color

RNA.celldict  %>% group_by(age, subclass_label_id) %>% dplyr::summarise(num=n()) %>% group_by(age) %>% dplyr::mutate(percentage=num/sum(num)) -> RNA.subclass.cellprop.age
RNA.celldict  %>% group_by(age_rep, subclass_label_id) %>% dplyr::summarise(num=n()) %>% group_by(age_rep) %>% dplyr::mutate(percentage=num/sum(num)) -> RNA.subclass.cellprop.age_rep

# Plot
pt_subclass.cellprop.age <- ggplot(RNA.subclass.cellprop.age, aes(x=age, y=percentage, fill=subclass_label_id)) + geom_area(alpha=0.6, size=0.5, colour="black")
pt_subclass.cellprop.age_theme <- pt_subclass.cellprop.age + scale_fill_manual(values=RNA_subclass_color$subclass_color) + scale_x_continuous(breaks=c(2, 6, 12, 18), limits=c(2, 18)) +
  xlab("Age (months)") + ylab("Percentage") + theme_bw() + 
  theme(legend.title=element_text(size=12), legend.text=element_text(size=10), axis.text.x=element_text(colour="black", size=16), 
        axis.text.y=element_text(colour="black", size=16), axis.title=element_text(colour="black", size=18))
ggsave(pt_subclass.cellprop.age_theme, file="./figures/fig2/cellsubclass_proportion_age.pdf", width=24, height=8)


# cell proportion differences
library("speckle")
library("limma")
# Run propeller testing for cell type proportion differences between the two groups
# subclass
# for continuous variable
props <- getTransformedProps(RNA$subclass_label_id, RNA$sampleID, transform="logit")
age <- rep(c(2,6,12,18), each=2) 
des.age <- model.matrix(~age)
des.age
# model on the Transformed proportions 
fit <- lmFit(props$TransformedProps,des.age)
fit <- eBayes(fit, robust=TRUE)
topTable(fit, n=Inf)
# model on the proportions 
fit.prop <- lmFit(props$Proportions,des.age)
fit.prop <- eBayes(fit.prop, robust=TRUE)
topTable(fit.prop, n=Inf)

#propeller_test_fit <- as.data.frame(fit)
#propeller_test_fit.prop <- as.data.frame(fit.prop)
subclass_propeller_test_fit <- topTable(fit, n=Inf)
subclass_propeller_test_fit.prop <- topTable(fit.prop, n=Inf)
write.csv(subclass_propeller_test_fit, file="./rds/RNA/diff/after_integra/subclass_propeller_test_fit.csv")
write.csv(subclass_propeller_test_fit.prop, file="./rds/RNA/diff/after_integra/subclass_propeller_test_fit.prop.csv")




