# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 3
# LSC peak fold change vs Bulk/Blast peak fold change

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")

lsc_cases<-unique(sample_info$patient[sample_info$type=="LSC"])
lsc_sample_info<-sample_info[sample_info$patient%in%lsc_cases&sample_info$type%in%c("LSC","Blast","Bulk"),]

genotype_groups<-genotyping_info[,c("ID","genomic_bin")]
lsc_sample_info$genotype_groups<-genotype_groups$genomic_bin[match(lsc_sample_info$patient,genotype_groups$ID)]


# Do LSC REL vs DX differential accessibility

# peaks
# read ATAC count matrix
count_matrix<-readRDS("inputs/bulk_count_matrix_raw.rds")
rownames(count_matrix)<-paste(count_matrix$Chr,count_matrix$Start,count_matrix$End,sep="_")
count_matrix<-count_matrix[,4:length(colnames(count_matrix))]

## paired DEseq analysis of stable samples
lsc_sample_info_sub<-lsc_sample_info[lsc_sample_info$genotype_groups=="Stable"&lsc_sample_info$type=="LSC",]
lsc_sample_info_sub<-lsc_sample_info_sub[!is.na(lsc_sample_info_sub$name),]
deseq_matrix<-count_matrix[,lsc_sample_info_sub$name]

data_reference<-lsc_sample_info_sub[,c("name","patient","timepoint")]
colnames(data_reference)<-c("sample","patient","timepoint")
dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~timepoint+patient))

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
deseq_Results<-as.data.frame(deseq_Results)
saveRDS(deseq_Results, "outputs/lsc_analysis/LSC_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")



# Do Blast REL vs DX differential accessibility

# peaks
# read ATAC count matrix
count_matrix<-readRDS("inputs/bulk_count_matrix_raw.rds")
rownames(count_matrix)<-paste(count_matrix$Chr,count_matrix$Start,count_matrix$End,sep="_")
count_matrix<-count_matrix[,4:length(colnames(count_matrix))]

## paired DEseq analysis of stable samples
lsc_sample_info_sub<-lsc_sample_info[lsc_sample_info$genotype_groups=="Stable"&lsc_sample_info$type=="Blast",]
lsc_sample_info_sub<-lsc_sample_info_sub[!is.na(lsc_sample_info_sub$name),]
deseq_matrix<-count_matrix[,lsc_sample_info_sub$name]

data_reference<-lsc_sample_info_sub[,c("name","patient","timepoint")]
colnames(data_reference)<-c("sample","patient","timepoint")
dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~timepoint+patient))

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
deseq_Results<-as.data.frame(deseq_Results)
saveRDS(deseq_Results, "outputs/lsc_analysis/Blast_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")



# Do Bulk REL vs DX differential accessibility

# peaks
# read ATAC count matrix
count_matrix<-readRDS("inputs/bulk_count_matrix_raw.rds")
rownames(count_matrix)<-paste(count_matrix$Chr,count_matrix$Start,count_matrix$End,sep="_")
count_matrix<-count_matrix[,4:length(colnames(count_matrix))]

## paired DEseq analysis of stable samples
lsc_sample_info_sub<-lsc_sample_info[lsc_sample_info$genotype_groups=="Stable"&lsc_sample_info$type=="Bulk",]
lsc_sample_info_sub<-lsc_sample_info_sub[!is.na(lsc_sample_info_sub$name),]
deseq_matrix<-count_matrix[,lsc_sample_info_sub$name]

data_reference<-lsc_sample_info_sub[,c("name","patient","timepoint")]
colnames(data_reference)<-c("sample","patient","timepoint")
dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~timepoint+patient))

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
deseq_Results<-as.data.frame(deseq_Results)
saveRDS(deseq_Results, "outputs/lsc_analysis/Bulk_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")



# plotting

lsc_diff_exp<-readRDS("outputs/lsc_analysis/LSC_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")
blast_diff_exp<-readRDS("outputs/lsc_analysis/Blast_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")
bulk_diff_exp<-readRDS("outputs/lsc_analysis/Bulk_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")

# LSC vs blast
regions<-intersect(rownames(lsc_diff_exp), rownames(blast_diff_exp))

plot_input<-data.frame(blast=blast_diff_exp$log2FoldChange[match(regions, rownames(blast_diff_exp))],
                       lsc=lsc_diff_exp$log2FoldChange[match(regions, rownames(lsc_diff_exp))])


plt<-ggplot(plot_input, aes(x=lsc, y=blast)) +
  ggrastr::rasterise(geom_point(size=0.1, alpha=1, color="grey40"), dpi=300) +
  theme_classic() +
  stat_cor(method = "pearson", label.sep = "\n", size=3) +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  xlab("LSC Paired Log Fold Change") +
  ylab("Non-LSC Paired Log Fold Change") +
  ggtitle("LSC vs Non-LSC Relapse vs Diagnosis Fold Change")

pdf("outputs/figure_3_LSC/LSC_vs_Blast_peak_fold_change_dotplot.pdf", width = 3, height = 4)
plt
dev.off()


# LSC vs bulk
regions<-intersect(rownames(lsc_diff_exp), rownames(bulk_diff_exp))

plot_input<-data.frame(bulk=bulk_diff_exp$log2FoldChange[match(regions, rownames(bulk_diff_exp))],
                       lsc=lsc_diff_exp$log2FoldChange[match(regions, rownames(lsc_diff_exp))])


plt<-ggplot(plot_input, aes(x=lsc, y=bulk)) +
  ggrastr::rasterise(geom_point(size=0.1, alpha=1, color="grey40"), dpi=300) +
  theme_classic() +
  stat_cor(method = "pearson", label.sep = "\n", size=3) +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  xlab("LSC Paired Log Fold Change") +
  ylab("Bulk Paired Log Fold Change") +
  ggtitle("LSC vs Bulk Relapse vs Diagnosis Fold Change")

pdf("outputs/figure_3_LSC/LSC_vs_Bulk_peak_fold_change_dotplot.pdf", width = 3, height = 4)
plt
dev.off()


# Blast vs bulk
regions<-intersect(rownames(bulk_diff_exp), rownames(blast_diff_exp))

plot_input<-data.frame(blast=blast_diff_exp$log2FoldChange[match(regions, rownames(blast_diff_exp))],
                       bulk=bulk_diff_exp$log2FoldChange[match(regions, rownames(bulk_diff_exp))])


plt<-ggplot(plot_input, aes(x=blast, y=bulk)) +
  ggrastr::rasterise(geom_point(size=0.1, alpha=1, color="grey40"), dpi=300) +
  theme_classic() +
  stat_cor(method = "pearson", label.sep = "\n", size=3) +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  xlab("Non-LSC Paired Log Fold Change") +
  ylab("Bulk Paired Log Fold Change") +
  ggtitle("Non-LSC vs Bulk Relapse vs Diagnosis Fold Change")

pdf("outputs/figure_3_LSC/Bulk_vs_Blast_peak_fold_change_dotplot.pdf", width = 3, height = 4)
plt
dev.off()
