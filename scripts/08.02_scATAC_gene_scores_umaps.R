# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Plot gene enrichments on umaps for each patient
#

library(ArchR)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")


# load and plot seaprate archr projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(ArchRProjects)<-unique(scATAC_ref$patient)

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
  
  ArchRProjects[[p]]<-addImputeWeights(ArchRProjects[[p]])
  
  all_genes<-getFeatures(ArchRProj = ArchRProjects[[p]], useMatrix = "GeneScoreMatrix")
  
  genes<-c("GATA1","GATA2","CEBPA","HOXA1","RUNX1","SPI1","MYC","MYB","FOS","JUN","IRF1","IRF2")
  
  genes_to_plot<-genes
  
  for(g in genes_to_plot){
    plt <- plotEmbedding(ArchRProj = ArchRProjects[[p]], plotAs = "points",
                         colorBy = "GeneScoreMatrix", 
                         name = g, 
                         #pal = ArchRPalettes$fireworks,
                         embedding = "UMAP_combined",
                         rastr = T, dpi=250,
                         xlim=xlims[[p]],ylim=ylims[[p]])
    
    #visualize two plots above side-by-side
    pdf(paste0("outputs/scATAC_plots/gene_scores_umaps/",p,"_",g,"_gene_score_umap.pdf"), width = 5, height = 5)
    print(plt)
    dev.off()
  }
}