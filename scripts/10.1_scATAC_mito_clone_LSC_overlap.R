# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 5
# Identify fraction of LSCs in each mito clone
#
#
# Workflow for this script is as follows:
# 1. 
# 2. 
# 3. 

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

samples<-c("SU142", "SU360", "SU892")


for(p in samples){
  lsc_classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))
  lsc_classifications$ID<-paste(lsc_classifications$timepoint, strsplit(rownames(lsc_classifications),"#") %>% sapply('[',2),sep="_")
  mito_clones<-readRDS(paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  
  mito_clones$LSC<-lsc_classifications$classification[match(mito_clones$ID,lsc_classifications$ID)]
  
  plt_input<-as.data.frame(table(mito_clones[,c("cluster","LSC")]))
  
  plt<-ggplot(plt_input, aes(fill=LSC, y=Freq, x=cluster)) + 
    geom_bar(position="fill", stat="identity", color="black", size=1) +
    scale_fill_manual(values = c("orange","yellow")) +
    theme_classic() +
    xlab("Mitochondrial Clone") +
    ylab("Cell Type Frequency") +
    ggtitle(p)
  print(plt)
}
