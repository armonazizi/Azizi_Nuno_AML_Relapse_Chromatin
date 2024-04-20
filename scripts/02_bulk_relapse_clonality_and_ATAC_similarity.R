# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 2
# Analysis of diagnosis vs relapse chromatin similarity across samples in different clonality groups

library(circlize)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(tibble)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")

source("scripts/00_global_ATAC_functions.R")

# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")
bulk_samples<-sample_info$ID[sample_info$is_bulk==1]
bulk_sample_info<-sample_info[sample_info$ID%in%bulk_samples,]
genotype_groups<-genotyping_info[,c("ID","genomic_bin")]
bulk_sample_info$genotype_groups<-genotype_groups$genomic_bin[match(bulk_sample_info$patient,genotype_groups$ID)]

# read ATAC count matrix
count_matrix<-readRDS("inputs/bulk_count_matrix_raw.rds")
rownames(count_matrix)<-paste(count_matrix$Chr,count_matrix$Start,count_matrix$End,sep="_")
count_matrix<-count_matrix[,4:length(colnames(count_matrix))]

# Normalize ATAC count matrix
count_matrix_normalized<-normalize_count_matrix(count_matrix)

# Read and normalize ATAC gene scores matrix
gene_scores<-readRDS("inputs/bulk_gene_scores.rds")
gene_scores<-column_to_rownames(gene_scores,"Gene")
gene_scores_normalized<-normalize_count_matrix(gene_scores)

count_matrix_normalized<-count_matrix_normalized[,bulk_sample_info$name]
gene_scores_normalized<-gene_scores_normalized[,bulk_sample_info$name]

## calculate dx-rel similarity

count_matrix_collapsed<-collapse_count_matrix(count_matrix_normalized, bulk_sample_info)
gene_scores_collapsed<-collapse_count_matrix(gene_scores_normalized, bulk_sample_info)

# Find ATAC similarity for each patient
correlations<-cor(count_matrix_collapsed+1)
corr_df<-data.frame()
for(p in unique(bulk_sample_info$patient)){
  corr<-correlations[bulk_sample_info$sample[bulk_sample_info$patient==p&bulk_sample_info$timepoint=="DX"][1],
                     bulk_sample_info$sample[bulk_sample_info$patient==p&bulk_sample_info$timepoint=="REL"][1]]
  corr_df[p,"Correlation"]<-corr
}
corr_df$genotype_group<-bulk_sample_info$genotype_groups[match(rownames(corr_df), bulk_sample_info$sample)]
corr_df<-corr_df[complete.cases(corr_df),]
corr_df$genotype_group[corr_df$genotype_group=="Gain_Loss"]<-"Gain+Loss"
corr_df$genotype_group<-factor(corr_df$genotype_group, levels = c("Stable","Gain","Loss","Gain+Loss"))

corr_df %>% group_by(genotype_group) %>% summarise(mean = mean(Correlation), n = n())

col = c(Stable = "#EE3A2D", Gain = "#2270B6", Loss="#70ACD4", `Gain+Loss`="#193D6C")

# do plotting
plt<-ggplot(corr_df, aes(x=genotype_group, y=Correlation, fill=genotype_group)) +
  geom_violin() +
  geom_jitter(width = 0.2) +
  scale_fill_manual(values=col) +
  theme_classic() +
  scale_y_continuous(limits=c(0.6,1)) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("DX vs REL Global Chromatin Similarity") +
  #ggtitle("DX vs REL peak accessibility similarity") +
  stat_compare_means(method="anova")
  #stat_compare_means(method="t.test", paired = FALSE, comparisons = list(c("Gain","Gain_Loss"),c("Gain","Stable")))

print(plt)

pdf("outputs/figure_2/clonal_groups_global_similarity.pdf", width = 4, height = 4)
print(plt)
dev.off()
