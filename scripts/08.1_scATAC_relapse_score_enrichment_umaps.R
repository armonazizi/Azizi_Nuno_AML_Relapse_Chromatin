# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Plot relapse signature on single cell umaps
#
# Workflow for this script is as follows:
# 1. score each cell based on relapse gene accessibility signature
# 2. Plot umaps with signature enrichment
# 3. 

library(ArchR)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")


# load and plot seaprate archr projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(ArchRProjects)<-unique(scATAC_ref$patient)

# do umaps
for(p in unique(scATAC_ref$patient)){
  message(p)
  
  relapse_signature<-readRDS("outputs/bulk_analysis/REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")
  
  relapse_signature<-relapse_signature[complete.cases(relapse_signature),]
  
  top_regions<-rownames(relapse_signature)[rev(order(relapse_signature$log2FoldChange))][1:500]
  bottom_regions<-rownames(relapse_signature)[order(relapse_signature$log2FoldChange)][1:500]
  
  top_regions<-GRanges(seqnames = strsplit(top_regions,"_") %>% sapply('[',1),
                       ranges=IRanges(start=as.numeric(strsplit(top_regions,"_") %>% sapply('[',2)), end=as.numeric(strsplit(top_regions,"_") %>% sapply('[',3))))
  bottom_regions<-GRanges(seqnames = strsplit(bottom_regions,"_") %>% sapply('[',1),
                          ranges=IRanges(start=as.numeric(strsplit(bottom_regions,"_") %>% sapply('[',2)), end=as.numeric(strsplit(bottom_regions,"_") %>% sapply('[',3))))
  
  ArchRProjects[[p]]<-addPeakSet(
    ArchRProj = ArchRProjects[[p]],
    peakSet = c(top_regions, bottom_regions),
    genomeAnnotation = getGenomeAnnotation(ArchRProjects[[p]]),
    force = T
  )
  
  ArchRProjects[[p]] <- addPeakMatrix(ArchRProjects[[p]])
  
  ArchRProjects[[p]]<-addFeatureCounts(
    ArchRProj = ArchRProjects[[p]],
    features = top_regions, 
    name = "rel_up")
  
  ArchRProjects[[p]]<-addFeatureCounts(
    ArchRProj = ArchRProjects[[p]],
    features = bottom_regions, 
    name = "rel_dn")
  
  ArchRProjects[[p]]@cellColData[,"rel_sig"]<-ArchRProjects[[p]]@cellColData[,"rel_upRatio"]-ArchRProjects[[p]]@cellColData[,"rel_dnRatio"] 
  
  
  #plot UMAP
  p1 <- plotEmbedding(ArchRProj = ArchRProjects[[p]], 
                      colorBy = "cellColData", 
                      name = "rel_sig", 
                      #pal = ArchRPalettes$fireworks,
                      embedding = "UMAP_postfilter",
                      rastr = T, dpi=250)
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_umap_relapse_signature_peaks_enrichment.pdf"), width = 5, height = 5)
  print(p1)
  dev.off()
}



### Do umaps using combined umap coordinates ###

# Load combined ArchR project
ArchR_comb<-loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/all_samples_combined"))

xlims<-list(SU142=c(-10,-3),
            SU360=c(-3,10),
            SU484=c(-2,14),
            SU892=c(-15,-3))

ylims<-list(SU142=c(3,10),
            SU360=c(-15,-2),
            SU484=c(-1,15),
            SU892=c(-9,3))

# do umaps
for(p in unique(scATAC_ref$patient)){
  message(p)
  
  ArchRProjects[[p]]@embeddings$UMAP_combined<-ArchRProjects[[p]]@embeddings$UMAP_postfilter
  ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df<-ArchR_comb@embeddings$UMAP_postfilter@listData$df[rownames(ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df),]
  
  ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df$`IterativeLSI_postfilter#UMAP_Dimension_1`[ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df$`IterativeLSI_postfilter#UMAP_Dimension_1`<xlims[[p]][1]]<-NA
  ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df$`IterativeLSI_postfilter#UMAP_Dimension_1`[ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df$`IterativeLSI_postfilter#UMAP_Dimension_1`>xlims[[p]][2]]<-NA
  ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df$`IterativeLSI_postfilter#UMAP_Dimension_2`[ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df$`IterativeLSI_postfilter#UMAP_Dimension_2`<ylims[[p]][1]]<-NA
  ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df$`IterativeLSI_postfilter#UMAP_Dimension_2`[ArchRProjects[[p]]@embeddings$UMAP_combined@listData$df$`IterativeLSI_postfilter#UMAP_Dimension_2`>ylims[[p]][2]]<-NA
  
  relapse_signature<-readRDS("outputs/bulk_analysis/REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")
  
  relapse_signature<-relapse_signature[complete.cases(relapse_signature),]
  
  top_regions<-rownames(relapse_signature)[rev(order(relapse_signature$log2FoldChange))][1:500]
  bottom_regions<-rownames(relapse_signature)[order(relapse_signature$log2FoldChange)][1:500]
  
  top_regions<-GRanges(seqnames = strsplit(top_regions,"_") %>% sapply('[',1),
                       ranges=IRanges(start=as.numeric(strsplit(top_regions,"_") %>% sapply('[',2)), end=as.numeric(strsplit(top_regions,"_") %>% sapply('[',3))))
  bottom_regions<-GRanges(seqnames = strsplit(bottom_regions,"_") %>% sapply('[',1),
                          ranges=IRanges(start=as.numeric(strsplit(bottom_regions,"_") %>% sapply('[',2)), end=as.numeric(strsplit(bottom_regions,"_") %>% sapply('[',3))))
  
  ArchRProjects[[p]]<-addPeakSet(
    ArchRProj = ArchRProjects[[p]],
    peakSet = c(top_regions, bottom_regions),
    genomeAnnotation = getGenomeAnnotation(ArchRProjects[[p]]),
    force = T
  )
  
  ArchRProjects[[p]] <- addPeakMatrix(ArchRProjects[[p]])
  
  ArchRProjects[[p]]<-addFeatureCounts(
    ArchRProj = ArchRProjects[[p]],
    features = top_regions, 
    name = "rel_up")
  
  ArchRProjects[[p]]<-addFeatureCounts(
    ArchRProj = ArchRProjects[[p]],
    features = bottom_regions, 
    name = "rel_dn")
  
  ArchRProjects[[p]]@cellColData[,"rel_sig"]<-ArchRProjects[[p]]@cellColData[,"rel_upRatio"]-ArchRProjects[[p]]@cellColData[,"rel_dnRatio"] 
  
  
  #plot UMAP
  p1 <- plotEmbedding(ArchRProj = ArchRProjects[[p]],
                      plotAs = "points",
                      colorBy = "cellColData", 
                      name = "rel_sig", 
                      #pal = ArchRPalettes$fireworks,
                      embedding = "UMAP_combined",
                      rastr = T, dpi=250,
                      xlim=xlims[[p]],ylim=ylims[[p]])
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_umap_relapse_signature_peaks_enrichment_combined_coords.pdf"), width = 5, height = 5)
  print(p1)
  dev.off()
  
  
  cells_to_plot<-rownames(ArchRProjects[[p]]@cellColData)[nchar(ArchRProjects[[p]]@cellColData$Sample)==5]
  
  #plot UMAP
  p1 <- plotEmbedding(ArchRProj = ArchRProjects[[p]], 
                      plotAs = "points",
                      highlightCells = cells_to_plot,
                      colorBy = "cellColData", 
                      name = "rel_sig", 
                      #pal = ArchRPalettes$fireworks,
                      embedding = "UMAP_combined",
                      rastr = T, dpi=250,
                      xlim=xlims[[p]],ylim=ylims[[p]])
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_umap_relapse_signature_peaks_enrichment_dx_highlight_combined_coords.pdf"), width = 5, height = 5)
  print(p1)
  dev.off()
}
 