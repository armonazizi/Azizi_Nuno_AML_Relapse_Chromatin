# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 3
# Bulk fold change vs LSC signature

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")


# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")

lsc_cases<-unique(sample_info$patient[sample_info$type=="LSC"])
lsc_sample_info<-sample_info[sample_info$patient%in%lsc_cases&sample_info$type%in%c("LSC","Blast","Bulk"),]

genotype_groups<-genotyping_info[,c("ID","genomic_bin")]
lsc_sample_info$genotype_groups<-genotype_groups$genomic_bin[match(lsc_sample_info$patient,genotype_groups$ID)]


# Do LSC vs non-LSC DX differential accessibility

# peaks
# read ATAC count matrix
count_matrix<-readRDS("inputs/bulk_count_matrix_raw.rds")
rownames(count_matrix)<-paste(count_matrix$Chr,count_matrix$Start,count_matrix$End,sep="_")
count_matrix<-count_matrix[,4:length(colnames(count_matrix))]

## paired DEseq analysis of stable samples
lsc_sample_info_sub<-lsc_sample_info[lsc_sample_info$timepoint=="DX"&lsc_sample_info$type%in%c("LSC","Blast")&lsc_sample_info$genotype_groups=="Stable",]
lsc_sample_info_sub<-lsc_sample_info_sub[!is.na(lsc_sample_info_sub$name),]
deseq_matrix<-count_matrix[,lsc_sample_info_sub$name]

data_reference<-lsc_sample_info_sub[,c("name","patient","type")]
colnames(data_reference)<-c("sample","patient","type")
dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~type+patient))

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("type","LSC","Blast"))
deseq_Results<-as.data.frame(deseq_Results)
saveRDS(deseq_Results, "outputs/lsc_analysis/DX_LSC_vs_Blast_peak_accessibility_deseq_res_stable_cases.rds")



# plotting

lsc_sig<-readRDS("outputs/lsc_analysis/DX_LSC_vs_Blast_peak_accessibility_deseq_res_stable_cases.rds")
bulk_diff_exp<-readRDS("outputs/lsc_analysis/Bulk_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")

# Bulk relapse accessibility vs LSC signature
regions<-intersect(rownames(lsc_sig), rownames(bulk_diff_exp))

plot_input<-data.frame(bulk=bulk_diff_exp$log2FoldChange[match(regions, rownames(bulk_diff_exp))],
                       lsc=lsc_sig$log2FoldChange[match(regions, rownames(lsc_sig))])


plt<-ggplot(plot_input, aes(x=bulk, y=lsc)) +
  ggrastr::rasterise(geom_point(size=0.1, alpha=1, color="grey40"), dpi=300) +
  theme_classic() +
  stat_cor(method = "pearson", label.sep = "\n", size=3) +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  xlab("Bulk Relapse Log Fold Change") +
  ylab("LSC Signature Log Fold Change") +
  ggtitle("Bulk Relapse Change vs LSC Signature")

pdf("outputs/figure_3_LSC/Bulk_peak_fold_change_vs_LSC_signature_dotplot.pdf", width = 3, height = 4)
plt
dev.off()