# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 5
# identify the fraction of cell states in each mito clone
#

library(ArchR)
library(openxlsx)
library(Seurat)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference_hg19.rds")
scATAC_ref$timepoint<-rep("DX", nrow(scATAC_ref))
scATAC_ref$timepoint[nchar(scATAC_ref$ID)==6]<-"REL"

for(p in unique(scATAC_ref$patient)){
  assignments<-readRDS(paste0("outputs/scATAC/AML_scATAC_cell_classifications/",p,"_KNN_classifications.rds"))
  assignments$timepoint<-rep("DX",nrow(assignments))
  assignments$timepoint[assignments$Sample=="reference"]<-"reference"
  assignments$timepoint[nchar(assignments$Sample)==6]<-"REL"
  assignments$ID<-paste(assignments$timepoint, strsplit(rownames(assignments),"#") %>% sapply('[',2),sep="_")
  
  mito_clones<-readRDS(paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  
  mito_clones$classification<-assignments$AML_classification[match(mito_clones$ID,assignments$ID)]
  
  mito_clones$clone<-paste(mito_clones$cluster,mito_clones$timepoint,sep="_")
  
  plt_input<-as.data.frame(table(mito_clones[,c("cluster","classification")]))
  
  celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono/DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
  
  plt<-ggplot(plt_input, aes(fill=classification, y=Freq, x=cluster)) + 
    geom_bar(position="fill", stat="identity", color="black", size=1) +
    scale_fill_manual(values = celltype_pal) +
    theme_classic() +
    xlab("Mitochondrial Clone") +
    ylab("Cell Type Frequency") +
    ggtitle(p)
  print(plt)
}
