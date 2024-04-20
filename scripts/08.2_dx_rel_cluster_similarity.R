# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# 
# Cluster DX and REL separately, then compare cluster epigenetic similarity
#

library(ArchR)
library(ggpubr)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load and plot seaprate archr projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(ArchRProjects)<-unique(scATAC_ref$patient)

for(p in unique(scATAC_ref$patient)){
  
  message(p)
  
  archr_dx<-loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",p,"_DX"))
  archr_rel<-loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",p,"_REL"))
  
  ArchRProjects[[p]]@cellColData$timepoint<-rep("DX",nrow(ArchRProjects[[p]]@cellColData))
  ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
  
  ArchRProjects[[p]]@cellColData$timepoint_clusters<-rep(NA, nrow(ArchRProjects[[p]]@cellColData))
  
  ArchRProjects[[p]]@cellColData$timepoint_clusters[ArchRProjects[[p]]@cellColData$timepoint=="DX"]<-paste0("DX_",archr_dx@cellColData$Clusters[match(rownames(ArchRProjects[[p]]@cellColData[ArchRProjects[[p]]@cellColData$timepoint=="DX",]), rownames(archr_dx@cellColData))])
  ArchRProjects[[p]]@cellColData$timepoint_clusters[ArchRProjects[[p]]@cellColData$timepoint=="REL"]<-paste0("REL_",archr_rel@cellColData$Clusters[match(rownames(ArchRProjects[[p]]@cellColData[ArchRProjects[[p]]@cellColData$timepoint=="REL",]), rownames(archr_rel@cellColData))])
  
  
  #plot UMAP
  p1 <- plotEmbedding(ArchRProj = archr_dx, 
                      colorBy = "cellColData", 
                      name = "Clusters", 
                      #pal = c("#272E6A","#D51F26"), 
                      embedding = "UMAP", 
                      size = 1,
                      rastr = T, dpi=250)
  
  p2 <- plotEmbedding(ArchRProj = archr_rel, 
                      colorBy = "cellColData", 
                      name = "Clusters", 
                      #pal = c("#272E6A","#D51F26"), 
                      embedding = "UMAP",
                      size = 1,
                      rastr = T, dpi=250)
  
  
  peak_mtx<-getMatrixFromProject(ArchRProj = ArchRProjects[[p]], useMatrix = "PeakMatrix")
  
  peak_mtx_collapsed<-as.data.frame(matrix(nrow=nrow(peak_mtx), ncol=0))
  
  for(c in unique(ArchRProjects[[p]]@cellColData$timepoint_clusters)){
    
    cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clusters==c]
    
    t<-strsplit(c,"_")[[1]][2]
    
    #cells_of_interest<-sample(cells_of_interest, 200, replace = T)
    
    pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
    colnames(pseudobulk)<-c(c)
    
    peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
  }
  
  peak_mtx_collapsed<-edgeR::cpm(peak_mtx_collapsed)
  
  sample_corr<-cor(peak_mtx_collapsed)
  
  annotations<-data.frame(timepoint=strsplit(colnames(sample_corr),"_") %>% sapply('[',1),
                          cluster=strsplit(colnames(sample_corr),"_") %>% sapply('[',2))
  cluster_colors<-ArchRPalettes$stallion[1:length(unique(annotations$cluster))]
  names(cluster_colors)<-sort(unique(annotations$cluster))
  timepoint_colors<-c("darkblue","darkred")
  names(timepoint_colors)<-c("DX","REL")
  anno_colors<-list(cluster=cluster_colors,
                    timepoint=timepoint_colors)
  
  p3<-Heatmap(sample_corr, 
              col=viridis(10),
              #col=colorRamp2(seq(0.5,1,length.out = 10), viridis(10)),
              show_row_names=FALSE, 
              show_column_names=T, 
              #show_column_dend = FALSE,
              top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"), 
              #left_annotation = rowAnnotation(df=annotations_mod, col=anno_colors, show_annotation_name=F),
              heatmap_legend_param = list(title = "Pearson"),
              column_title=paste(p,"Cluster Epigenetic Similarity"),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8))
  
  ht_grob <- grid.grabExpr(draw(p3)) 
  
  
  pdf(paste0("outputs/scATAC_plots/cluster_similarity/",p,"_dx_rel_cluster_similarity.pdf"), width = 12, height = 4)
  print(ggarrange(p1, p2, ht_grob, nrow=1,widths = c(1,1,1.3)))
  dev.off()
  
}
