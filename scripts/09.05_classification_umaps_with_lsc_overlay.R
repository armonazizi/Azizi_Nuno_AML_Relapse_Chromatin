# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Umap of projected celltypes with LSC overlay

library(ArchR)
library(openxlsx)
library(Seurat)
library(ggpubr)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")

ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(ArchRProjects)<-unique(scATAC_ref$patient)



for(p in unique(scATAC_ref$patient)){
  
  # label LSCs
  if(file.exists(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))){
    lsc_classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))
    
    ArchRProjects[[p]]@cellColData$CellType<-lsc_classifications$classification[match(rownames(ArchRProjects[[p]]@cellColData),rownames(lsc_classifications))]
  }else{
    next
  }
  
  # Label projected celltypes
  classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  ArchRProjects[[p]]@cellColData$projection_classification<-classifications$AML_classification[match(rownames(ArchRProjects[[p]]@cellColData),rownames(classifications))]
  
  plotEmbedding(ArchRProj = ArchRProjects[[p]], 
                colorBy = "cellColData", 
                name = "projection_classification",
                embedding = "UMAP_postfilter", size=0.5,
                labelMeans=T)
  
  coldata<-as.data.frame(ArchRProjects[[p]]@cellColData)
  
  plt_input<-as.data.frame(ArchRProjects[[p]]@embeddings$UMAP_postfilter$df)
  colnames(plt_input)<-c("UMAP_1","UMAP_2")
  plt_input$LSC<-coldata$CellType[match(rownames(plt_input), rownames(plt_input))]
  plt_input$classification<-coldata$projection_classification[match(rownames(plt_input), rownames(plt_input))]
  
  celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono-like","DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
  
  ggPoint(plt_input$UMAP_1,plt_input$UMAP_2, color = plt_input$classification, pal = celltype_pal, )
  
  plt<-ggplot(plt_input, aes(x=UMAP_1,y=UMAP_2,color=classification)) +
    ggrastr::rasterise(geom_point(size=0.5), dpi=300) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = celltype_pal) +
    ggrastr::rasterise(geom_point(data=plt_input[plt_input$LSC=="LSC",],aes(x=UMAP_1,y=UMAP_2), color="red",size=0.5), dpi=300) +
    ggtitle(p)
  
  pdf(paste0("outputs/scATAC_plots/",p,"_projection_classification_dotplot_LSC_overlay.pdf"), width = 5, height = 3.5, useDingbats = FALSE)
  print(plt)
  dev.off()
  
}


# Map LSCs onto projection umap

for(p in unique(scATAC_ref$patient)){
  
  assignments<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  
  # label LSCs
  if(file.exists(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))){
    lsc_classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))
    assignments$CellType<-lsc_classifications$classification[match(rownames(assignments),rownames(lsc_classifications))]
    assignments$CellType[is.na(assignments$CellType)]<-"reference"
  }else{
    next
  }
  
  assignments$CellType<-factor(assignments$CellType, levels = c("reference","Blast","LSC"))
  assignments<-assignments[order(assignments$CellType),]
  
  pal<-c("lightgrey","darkgreen","darkorange")
  names(pal)<-c("reference","LSC","Blast")
  
  plt<-ggplot(assignments, aes(x=X1,y=X2,color=CellType)) +
    ggrastr::rasterise(geom_point(size=0.5), dpi=300) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = pal) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(p)
  
  pdf(paste0("outputs/scATAC_plots/",p,"_LSI_projection_umap_LSC_color.pdf"), width = 5, height = 4, useDingbats = FALSE)
  print(plt)
  dev.off()
  
}
