# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Celltype projection fractions within LSC compartment vs non-LSC compartment

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

samples<-c("SU142", "SU360", "SU892")

for(p in samples){
  lsc_classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))
  lsc_classifications$ID<-rownames(lsc_classifications)
  assignments<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  assignments$ID<-rownames(assignments)
  assignments<-assignments[assignments$AML_classification!="reference",]
  
  assignments$LSC<-lsc_classifications$classification[match(assignments$ID,lsc_classifications$ID)]
  
  plt_input<-as.data.frame(table(assignments[,c("AML_classification","LSC")]))
  
  celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono-like","DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
  
  plt<-ggplot(plt_input, aes(fill=AML_classification, y=Freq, x=LSC)) + 
    geom_bar(position="fill", stat="identity", color="black", size=1) +
    scale_fill_manual(values = celltype_pal) +
    theme_classic() +
    xlab("LSC Classification") +
    ylab("Cell Type Frequency") +
    ggtitle(p)
  print(plt)
}


for(p in samples){
  lsc_classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))
  lsc_classifications$ID<-rownames(lsc_classifications)
  assignments<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  assignments$ID<-rownames(assignments)
  assignments<-assignments[assignments$AML_classification!="reference",]
  
  assignments$LSC<-lsc_classifications$classification[match(assignments$ID,lsc_classifications$ID)]
  
  assignments$timepoint<-rep("DX", nrow(assignments))
  assignments$timepoint[nchar(assignments$Sample)==6]<-"REL"
  
  assignments$LSC_timepoint<-paste(assignments$timepoint, assignments$LSC, sep=" ")
  
  plt_input<-as.data.frame(table(assignments[,c("AML_classification","LSC_timepoint")]))
  
  celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono-like","DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
  
  plt<-ggplot(plt_input, aes(fill=AML_classification, y=Freq, x=LSC_timepoint)) + 
    geom_bar(position="fill", stat="identity", color="black", size=0) +
    scale_fill_manual(values = celltype_pal) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("LSC Classification") +
    ylab("Cell Type Frequency") +
    ggtitle(p)
  
  pdf(paste0("outputs/scATAC_plots/",p,"_LSC_vs_projection_classifications_barplot.pdf"), width = 5, height = 3.5, useDingbats = FALSE)
  print(plt)
  dev.off()
  
}
