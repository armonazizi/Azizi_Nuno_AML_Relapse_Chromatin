# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 5
# Differential motif accessibility between mito subclones
#
# 1. calculate marker peaks for subclones
# 2. identify motifs enriched in each subclone

library(ArchR)
library(openxlsx)
library(Seurat)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
scATAC_ref$timepoint<-rep("DX", nrow(scATAC_ref))
scATAC_ref$timepoint[nchar(scATAC_ref$ID)==6]<-"REL"

ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(ArchRProjects)<-unique(scATAC_ref$patient)


for(p in unique(scATAC_ref$patient)){
  
  message(p)
  
  
  # label classification celltypes
  mito_clones<-readRDS(paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  
  mito_clones$ID[mito_clones$timepoint=="DX"]<-paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="DX"], "#",
                                                      sapply(strsplit(mito_clones$ID[mito_clones$timepoint=="DX"],"_"),'[',2))
  
  mito_clones$ID[mito_clones$timepoint=="REL"]<-paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="REL"], "#",
                                                       sapply(strsplit(mito_clones$ID[mito_clones$timepoint=="REL"],"_"),'[',2))
  
  ArchRProjects[[p]]@cellColData$mito_clone<-as.character(mito_clones$cluster[match(rownames(ArchRProjects[[p]]@cellColData),mito_clones$ID)])
  ArchRProjects[[p]]@cellColData$mito_clone[is.na(ArchRProjects[[p]]@cellColData$mito_clone)]<-"Unk"
  
  ArchRProjects[[p]]@cellColData$timepoint<-rep("DX", nrow(ArchRProjects[[p]]@cellColData))
  ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
  
  ArchRProjects[[p]]@cellColData$timepoint_clone<-paste(ArchRProjects[[p]]@cellColData$mito_clone,
                                                        ArchRProjects[[p]]@cellColData$timepoint, sep="_")
  
  
  ArchRProjects[[p]]<-addPeakSet(
    ArchRProj = ArchRProjects[[p]],
    peakSet = readRDS(paste0("inputs/scATAC/scATAC_peaksets/",p,"_ArchR_peakset.rds")),
    genomeAnnotation = getGenomeAnnotation(ArchRProjects[[p]]),
    force = T
  )
  
  ArchRProjects[[p]]<-addPeakMatrix(ArchRProjects[[p]])
  
  ArchRProjects[[p]] <- addBgdPeaks(ArchRProjects[[p]])
  
  ArchRProjects[[p]] <- addMotifAnnotations(ArchRProj = ArchRProjects[[p]], motifSet = "cisbp", name = "Motif", force = T)
  
  markersPeaks <- getMarkerFeatures(
    ArchRProj = ArchRProjects[[p]], 
    useMatrix = "PeakMatrix", 
    groupBy = "timepoint_clone",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  heatmapPeaks <- markerHeatmap(
    seMarker = markersPeaks, 
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    transpose = TRUE
  )
  
  draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  
  
  enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = ArchRProjects[[p]],
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
  
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  
  
  
  
}