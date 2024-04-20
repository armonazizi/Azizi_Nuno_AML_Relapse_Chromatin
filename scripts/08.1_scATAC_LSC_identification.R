# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# ID LSCs in single cell data using bulk signatures
#
# Workflow for this script is as follows:
# 1. Using bulk differential signatures, generate LSC-specific peakset
# 2. Score single cells using the LSC peakset
# 3. Classify single cells
# 4. Do plotting and differential analysis
#
#
# May be able to use addFeatureCounts() from ArchR to count reads in features for scoring

library(ArchR)
library(openxlsx)
library(Seurat)
library(ggpubr)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")

AML_ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(AML_ArchRProjects)<-unique(scATAC_ref$patient)

# Add LSC scores for each patient
for(p in c("SU142", "SU360", "SU892")){
  message(p)
  
  for(tp in c("DX","REL")){
    lsc_deseq_res<-readRDS(paste0("../LSC_relapse_analysis/per_patient_LSC_vs_blast/",p,"_",tp,"_LSC_vs_Blast_peaks_differential_accessibility.rds"))
    
    lsc_deseq_res<-lsc_deseq_res[complete.cases(lsc_deseq_res),]
    
    top_regions<-rownames(lsc_deseq_res)[rev(order(lsc_deseq_res$log2FoldChange))][1:50]
    bottom_regions<-rownames(lsc_deseq_res)[order(lsc_deseq_res$log2FoldChange)][1:50]
    
    top_regions<-GRanges(seqnames = strsplit(top_regions,"_") %>% sapply('[',1),
                         ranges=IRanges(start=as.numeric(strsplit(top_regions,"_") %>% sapply('[',2)), end=as.numeric(strsplit(top_regions,"_") %>% sapply('[',3))))
    bottom_regions<-GRanges(seqnames = strsplit(bottom_regions,"_") %>% sapply('[',1),
                            ranges=IRanges(start=as.numeric(strsplit(bottom_regions,"_") %>% sapply('[',2)), end=as.numeric(strsplit(bottom_regions,"_") %>% sapply('[',3))))
    
    AML_ArchRProjects[[p]]<-addPeakSet(
      ArchRProj = AML_ArchRProjects[[p]],
      peakSet = c(top_regions, bottom_regions),
      genomeAnnotation = getGenomeAnnotation(AML_ArchRProjects[[p]]),
      force = T
    )
    
    AML_ArchRProjects[[p]] <- addPeakMatrix(AML_ArchRProjects[[p]])
    
    AML_ArchRProjects[[p]]<-addFeatureCounts(
      ArchRProj = AML_ArchRProjects[[p]],
      features = top_regions, 
      name = paste0("lsc_up_",tp))
    
    AML_ArchRProjects[[p]]<-addFeatureCounts(
      ArchRProj = AML_ArchRProjects[[p]],
      features = bottom_regions, 
      name = paste0("lsc_dn_",tp))
    
    AML_ArchRProjects[[p]]@cellColData[,paste0("lsc_score_",tp)]<-AML_ArchRProjects[[p]]@cellColData[,paste0("lsc_up_",tp,"Ratio")]-AML_ArchRProjects[[p]]@cellColData[,paste0("lsc_dn_",tp,"Ratio")] 
  }
  
  p1<-plotEmbedding(ArchRProj = AML_ArchRProjects[[p]], embedding = "UMAP_postfilter", colorBy = "cellColData",name="Sample", pal = c("#D51F26", "#272E6A"))
  p2<-plotEmbedding(ArchRProj = AML_ArchRProjects[[p]], embedding = "UMAP_postfilter", colorBy = "cellColData",name="lsc_up_DXRatio")
  p3<-plotEmbedding(ArchRProj = AML_ArchRProjects[[p]], embedding = "UMAP_postfilter", colorBy = "cellColData",name="lsc_up_RELRatio")

  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_umap_LSC_scores.pdf"), width = 12, height = 5)
  print(ggarrange(p1,p2,p3, nrow=1))
  dev.off()
}


# classify LSCs and save cell barcodes

LSC_fractions<-data.frame(SU142=c(0.58,0.58),
                          SU360=c(0.08,0.04),
                          SU892=c(0.1,0.1))
rownames(LSC_fractions)<-c("DX","REL")

for(p in c("SU142", "SU360", "SU892")){
  message(p)
  
  cell_data<-as.data.frame(AML_ArchRProjects[[p]]@cellColData)
  cell_data$timepoint<-rep("DX",nrow(cell_data))
  cell_data$timepoint[nchar(cell_data$Sample)==6]<-"REL"
  
  num_dx_lscs<-sum(cell_data$timepoint=="DX")*LSC_fractions["DX",p]
  cell_data_sub<-cell_data[cell_data$timepoint=="DX",]
  dx_barcodes<-rownames(cell_data_sub[rev(order(cell_data_sub[,"lsc_up_DXRatio"])),])[1:num_dx_lscs]
  
  num_rel_lscs<-sum(cell_data$timepoint=="REL")*LSC_fractions["REL",p]
  cell_data_sub<-cell_data[cell_data$timepoint=="REL",]
  rel_barcodes<-rownames(cell_data_sub[rev(order(cell_data_sub[,"lsc_up_RELRatio"])),])[1:num_rel_lscs]
  
  cell_data<-cell_data[c("Sample","timepoint","lsc_up_DXRatio","lsc_up_RELRatio")]
  cell_data$classification<-rep("Blast",nrow(cell_data))
  cell_data$classification[rownames(cell_data)%in%c(dx_barcodes,rel_barcodes)]<-"LSC"

  saveRDS(cell_data, paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))
}
