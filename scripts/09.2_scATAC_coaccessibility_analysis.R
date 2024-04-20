# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Coaccessibility tracks at DX and REL across select loci
# Identification of relapse specific co-accessibility
#


library(ArchR)
library(openxlsx)
library(Seurat)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
scATAC_ref$timepoint<-rep("DX", nrow(scATAC_ref))
scATAC_ref$timepoint[nchar(scATAC_ref$ID)==6]<-"REL"

AML_ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(AML_ArchRProjects)<-unique(scATAC_ref$patient)


p<-"SU484"

ArchRProjects[[p]]@cellColData$timepoint<-rep("DX", nrow(ArchRProjects[[p]]@cellColData))
ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"

ArchRProjects[[p]] <- addCoAccessibility(
  ArchRProj = ArchRProjects[[p]], 
  #cellsToUse = rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint=="DX"],
  reducedDims = "IterativeLSI"
)

dx_coaccessibility<-getCoAccessibility(ArchRProjects[[p]])

ArchRProjects[[p]] <- addCoAccessibility(
  ArchRProj = ArchRProjects[[p]], 
  cellsToUse = rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint=="REL"],
  reducedDims = "IterativeLSI"
)

rel_coaccessibility<-getCoAccessibility(ArchRProjects[[p]])


# co-accessibility stats

length(dx_coaccessibility$CoAccessibility)
length(rel_coaccessibility$CoAccessibility)

length(intersect(paste(dx_coaccessibility$CoAccessibility@seqnames,dx_coaccessibility$CoAccessibility@ranges@start,dx_coaccessibility$CoAccessibility@ranges@width,sep="_"),
                 paste(rel_coaccessibility$CoAccessibility@seqnames,rel_coaccessibility$CoAccessibility@ranges@start,rel_coaccessibility$CoAccessibility@ranges@width,sep="_")))



HOXA_locus<-GRanges(seqnames = c("chr7"), ranges = IRanges(start=c(27080000), end=c(27220000)))
HOXB_locus<-GRanges(seqnames = c("chr17"), ranges = IRanges(start=c(48520000), end=c(48640000)))
HIST_locus<-GRanges(seqnames = c("chr6"), ranges = IRanges(start=c(26000000), end=c(26300000)))
CDK6_locus<-GRanges(seqnames = c("chr7"), ranges = IRanges(start=c(92489398), end=c(93088969)))

gene<-"FLT3"
width<-100000
locus<-HOXA_locus

#gene
if(T) {
  plt1 <- plotBrowserTrack(
    ArchRProj = ArchRProjects[[p]], 
    groupBy = "timepoint", 
    pal = c(DX="darkblue", REL="darkred"),
    #region = HIST_locus,
    geneSymbol = gene,
    #useGroups = c("00_DX","00_REL","01_DX","01_REL","02_DX","02_REL","03_DX","03_REL","04_DX","04_REL","05_DX","05_REL"),
    #plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    #loops = dx_coaccessibility,
    upstream = width,
    downstream = width
  )
  
  plt2 <- plotBrowserTrack(
    ArchRProj = ArchRProjects[[p]], 
    groupBy = "timepoint",
    pal = c("darkblue", REL="darkred"),
    #region = HIST_locus,
    geneSymbol = gene,
    #useGroups = c("00_DX","00_REL","01_DX","01_REL","02_DX","02_REL","03_DX","03_REL","04_DX","04_REL","05_DX","05_REL"),
    #plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    #loops = rel_coaccessibility,
    upstream = width,
    downstream = width, 
  )
  
  grid::grid.newpage()
  grid::grid.draw(plt1[[gene]])
  
  grid::grid.newpage()
  grid::grid.draw(plt2[[gene]]) 
}

#locus
if(T) {
  plt1 <- plotBrowserTrack(
    ArchRProj = ArchRProjects[[p]], 
    groupBy = "timepoint", 
    pal = c(DX="darkblue", REL="darkred"),
    region = locus,
    #geneSymbol = gene,
    #useGroups = c("00_DX","00_REL","01_DX","01_REL","02_DX","02_REL","03_DX","03_REL","04_DX","04_REL","05_DX","05_REL"),
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    loops = dx_coaccessibility,
    upstream = width,
    downstream = width, 
    title = paste(p, locus)
  )
  
  plt2 <- plotBrowserTrack(
    ArchRProj = ArchRProjects[[p]], 
    groupBy = "timepoint",
    pal = c("darkblue", REL="darkred"),
    region = locus,
    #geneSymbol = gene,
    #useGroups = c("00_DX","00_REL","01_DX","01_REL","02_DX","02_REL","03_DX","03_REL","04_DX","04_REL","05_DX","05_REL"),
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    loops = rel_coaccessibility,
    upstream = width,
    downstream = width, 
    title = paste(p, locus)
  )
  
  
  grid::grid.newpage()
  grid::grid.draw(plt1)
  
  grid::grid.newpage()
  grid::grid.draw(plt2) 
}

pdf("outputs/scATAC_plots/accessibility_tracks/SU360_HOXA_DX.pdf", width = 8, height=8, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(plt1)
dev.off()

pdf("outputs/scATAC_plots/accessibility_tracks/SU360_HOXA_REL.pdf", width = 8, height=8, useDingbats = F)
grid::grid.newpage()
grid::grid.draw(plt2)
dev.off()




plt1 <- plotBrowserTrack(
  ArchRProj = ArchRProjects[[p]], 
  groupBy = "timepoint", 
  pal = c(DX="darkblue", REL="darkred"),
  region = locus,
  #geneSymbol = gene,
  #useGroups = c("00_DX","00_REL","01_DX","01_REL","02_DX","02_REL","03_DX","03_REL","04_DX","04_REL","05_DX","05_REL"),
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  loops = getCoAccessibility(ArchRProjects[[p]]),
  upstream = width,
  downstream = width, 
  title = paste(p, locus)
)


#ggarrange(plt1[[gene]],plt2[[gene]], ncol=1)