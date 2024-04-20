# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 2
# kaplan meier analysis of time to relapse in stabford cohort


library(circlize)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(tibble)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")

source("scripts/00_global_ATAC_functions.R")

# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/Supplementary_PatientCharacteristics.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")

ttr_data<-genotyping_info[,c("ID","genomic_bin")]
ttr_data$TTR<-sample_info$Days.Between.Sample.Acquisitions[match(ttr_data$ID,sample_info$Patient)]
ttr_data$timepoint<-sample_info$Disease.Stage[match(ttr_data$ID,sample_info$Patient)]
ttr_data<-ttr_data[complete.cases(ttr_data),]
ttr_data<-ttr_data[ttr_data$timepoint=="Diagnosis",]

ttr_data$clonality<-ttr_data$genomic_bin
ttr_data$clonality[ttr_data$clonality!="Stable"]<-"Unstable Clonality"
ttr_data$clonality[ttr_data$clonality=="Stable"]<-"Stable Clonality"

fit<-survfit(Surv(TTR) ~ clonality, data = ttr_data)

plt<-ggsurvplot(fit, data = ttr_data,
                legend.labs = c("Stable Clonality", "Unstable Clonality"),
                pval = TRUE,
                palette = c("#EE3A2D", "#193D6C"), #### changed colors
                ylab="Relapse Probability",
                break.time.by=300, ### changed breaks
                xlab="Time (days)",
                xlim = c(0, 1500) ###added xlim
)

print(plt)

pdf("outputs/figure_2/stable_vs_nonstable_RFS_KM.pdf", width = 5, height = 4, )
print(plt[[1]])
dev.off()
