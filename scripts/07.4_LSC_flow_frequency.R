# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 3
# LSC frequency at DX and REL


setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/genotyping_final.xlsx")

lsc_freq<-sample_info[sample_info$sample!="SU654C",c("patient","timepoint","lsc_percent")]
lsc_freq<-lsc_freq[!duplicated(lsc_freq),]

genotype_groups<-genotyping_info[,c("Sample.ID","Genomic.Bin")]
lsc_freq$genotype_groups<-genotype_groups$Genomic.Bin[match(lsc_freq$patient,genotype_groups$Sample.ID)]
lsc_freq$genotype_groups[lsc_freq$genotype_groups!="Stable"]<-"Unstable"

lsc_freq<-lsc_freq[complete.cases(lsc_freq),]

plt<-ggplot(lsc_freq,aes(x=timepoint, y=lsc_percent,group=patient, color=genotype_groups)) +
  #geom_boxplot(mapping = aes(x=timepoint, y=lsc_percent), inherit.aes = F, size=1) +
  geom_line(size=0.5) +
  geom_point(shape=21, fill="white") +
  scale_color_manual(name = "Genetic\nEvolution", values=c("red","darkblue")) +
  xlab("") +
  ylab("LSC Frequency") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  stat_compare_means(method="t.test", paired = TRUE, comparisons = list(c("DX","REL")))

print(plt)

pdf("outputs/figure_3_LSC/lsc_flow_frequency_dotplot.pdf", width = 3, height = 3.5)
print(plt)
dev.off()
