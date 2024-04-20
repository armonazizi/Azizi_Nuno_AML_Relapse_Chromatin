# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 5
# Differential accessibility between mito subclones
#
# 1. do global chromatin similarity across subclones at DX and REL for each sample
# 2. compare specific predefined subclones for each sample to find differentially accessible features

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



## Global similarity across subclones
# do differential accessibility for each subpopulation - need to control for cell number
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
  
  
  
  #ArchRProjects[[p]]<-addGroupCoverages(ArchRProjects[[p]], "timepoint_clone")
  
  
  peak_mtx<-getMatrixFromProject(ArchRProj = ArchRProjects[[p]], useMatrix = "PeakMatrix")
  
  peak_mtx_collapsed<-as.data.frame(matrix(nrow=nrow(peak_mtx), ncol=0))
  
  for(c in unique(ArchRProjects[[p]]@cellColData$timepoint_clone)){
    if(c %in% c("Unk_DX","Unk_REL")){next}
    
    cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clone==c]
    
    t<-strsplit(c,"_")[[1]][2]
    
    #if(length(cells_of_interest)<70){next}
    if(length(cells_of_interest)/sum(ArchRProjects[[p]]@cellColData$timepoint==t)<0.05){next}
    
    cells_of_interest<-sample(cells_of_interest, 200, replace = T)
    
    pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
    colnames(pseudobulk)<-c(c)
    
    peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
  }
  
  if(p=="SU484"){
    c<-"02_DX"
    cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clone==c]
    pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
    colnames(pseudobulk)<-c(c)
    peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
  }
  
  peak_mtx_collapsed<-edgeR::cpm(peak_mtx_collapsed)
  
  sample_corr<-cor(peak_mtx_collapsed)
  
  annotations<-data.frame(clone=strsplit(colnames(sample_corr),"_") %>% sapply('[',1),
                          timepoint=strsplit(colnames(sample_corr),"_") %>% sapply('[',2))
  clone_colors<-ArchRPalettes$stallion[1:length(unique(annotations$clone))]
  names(clone_colors)<-sort(unique(annotations$clone))
  timepoint_colors<-c("darkblue","darkred")
  names(timepoint_colors)<-c("DX","REL")
  anno_colors<-list(clone=clone_colors,
                    timepoint=timepoint_colors)
  
  ht<-Heatmap(sample_corr, 
          col=viridis(10),
          #col=colorRamp2(seq(0.5,1,length.out = 10), viridis(10)),
          show_row_names=FALSE, 
          show_column_names=T, 
          #show_column_dend = FALSE,
          top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"), 
          #left_annotation = rowAnnotation(df=annotations_mod, col=anno_colors, show_annotation_name=F),
          heatmap_legend_param = list(title = "Pearson"),
          column_title=paste(p,"Subclone Epigenetic Similarity"),
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8))
  
  
  pdf(paste0("outputs/scATAC_plots/mito_subclone_accessibility/",p,"_mito_subclone_similarity.pdf"), width = 5, height = 4.5, useDingbats = FALSE)
  print(ht)
  dev.off()
  
  
  # Compare diagnosis clone similarity to relapse clone similarity
  
  plt_input<-reshape2::melt(sample_corr)
  plt_input$timepoint<-rep(NA, nrow(sample_corr))
  plt_input$timepoint[grepl("DX",plt_input$Var1)&grepl("DX",plt_input$Var2)]<-"DX"
  plt_input$timepoint[grepl("REL",plt_input$Var1)&grepl("REL",plt_input$Var2)]<-"REL"
  plt_input$timepoint[plt_input$Var1==plt_input$Var2]<-NA
  plt_input<-plt_input[complete.cases(plt_input),]
  
  p1 <- ggplot(plt_input, aes(x=timepoint, y=value, color=timepoint, fill=timepoint)) +
    geom_boxplot(color="black",outlier.alpha = 0, size=0.5) +
    ggrastr::rasterize(geom_jitter(width = 0.1, size=0.5, show.legend = F, color="black"), dpi=300) +
    theme_classic() +
    scale_color_manual(values=c("darkblue","darkred")) +
    scale_fill_manual(values=c("darkblue","darkred")) +
    theme(axis.text.x = element_text(color="black",size=10)) +
    xlab("") +
    ylab("Between-Clone Epigenetic Similarity") +
    stat_compare_means(method="t.test", comparisons = list(c("DX","REL")))
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/mito_subclone_accessibility/",p,"_DX_vs_REL_mito_subclone_similarity.pdf"), width = 3, height = 4)
  print(p1)
  dev.off()
  
}



# Shuffle labels, do similarity for each subpopulation 
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
  
  
  # SHUFFLE LABELS
  ArchRProjects[[p]]@cellColData$mito_clone<-sample(ArchRProjects[[p]]@cellColData$mito_clone, replace = F)
  
  ArchRProjects[[p]]@cellColData$timepoint<-rep("DX", nrow(ArchRProjects[[p]]@cellColData))
  ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
  
  ArchRProjects[[p]]@cellColData$timepoint_clone<-paste(ArchRProjects[[p]]@cellColData$mito_clone,
                                                        ArchRProjects[[p]]@cellColData$timepoint, sep="_")
  
  
  print(table(ArchRProjects[[p]]@cellColData$timepoint_clone))
  
  #ArchRProjects[[p]]<-addGroupCoverages(ArchRProjects[[p]], "timepoint_clone")
  
  
  peak_mtx<-getMatrixFromProject(ArchRProj = ArchRProjects[[p]], useMatrix = "PeakMatrix")
  
  peak_mtx_collapsed<-as.data.frame(matrix(nrow=nrow(peak_mtx), ncol=0))
  
  for(c in unique(ArchRProjects[[p]]@cellColData$timepoint_clone)){
    if(c %in% c("Unk_DX","Unk_REL")){next}
    
    cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clone==c]
    
    if(length(cells_of_interest)<50){next}
    
    cells_of_interest<-sample(cells_of_interest, 200, replace = T)
    
    pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
    colnames(pseudobulk)<-c(c)
    
    peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
  }
  
  peak_mtx_collapsed<-edgeR::cpm(peak_mtx_collapsed)
  
  sample_corr<-cor(peak_mtx_collapsed)
  
  annotations<-data.frame(clone=strsplit(colnames(sample_corr),"_") %>% sapply('[',1),
                          timepoint=strsplit(colnames(sample_corr),"_") %>% sapply('[',2))
  clone_colors<-ArchRPalettes$stallion[1:length(unique(annotations$clone))]
  names(clone_colors)<-sort(unique(annotations$clone))
  timepoint_colors<-c("darkblue","darkred")
  names(timepoint_colors)<-c("DX","REL")
  anno_colors<-list(clone=clone_colors,
                    timepoint=timepoint_colors)
  
  ht<-Heatmap(sample_corr, 
              col=viridis(10),
              #col=colorRamp2(seq(0.5,0.9,length.out = 10), viridis(10)),
              show_row_names=FALSE, 
              show_column_names=T, 
              #show_column_dend = FALSE,
              top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"), 
              #left_annotation = rowAnnotation(df=annotations_mod, col=anno_colors, show_annotation_name=F),
              heatmap_legend_param = list(title = "Pearson"),
              column_title=paste(p,"Subclone Epigenetic Similarity"),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8))
  
  
  pdf(paste0("outputs/scATAC_plots/mito_subclone_accessibility/",p,"_mito_subclone_similarity_labels_shuffled.pdf"), width = 5, height = 4.5, useDingbats = FALSE)
  print(ht)
  dev.off()
}




# Shuffle labels and build similarity distribution
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
  
  
  
  peak_mtx<-getMatrixFromProject(ArchRProj = ArchRProjects[[p]], useMatrix = "PeakMatrix")
  
  cor_set<-list()
  
  for(i in 1:100){
    if(i%%100==1){message(i)}
    
    ArchRProjects[[p]]@cellColData$timepoint<-rep("DX", nrow(ArchRProjects[[p]]@cellColData))
    ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
    
    ArchRProjects[[p]]@cellColData$timepoint_clone<-paste(ArchRProjects[[p]]@cellColData$mito_clone,
                                                          ArchRProjects[[p]]@cellColData$timepoint, sep="_")
    
    # SHUFFLE LABELS
    ArchRProjects[[p]]@cellColData$timepoint_clone[ArchRProjects[[p]]@cellColData$timepoint=="DX"]<-sample(ArchRProjects[[p]]@cellColData$timepoint_clone[ArchRProjects[[p]]@cellColData$timepoint=="DX"], replace = F)
    ArchRProjects[[p]]@cellColData$timepoint_clone[ArchRProjects[[p]]@cellColData$timepoint=="REL"]<-sample(ArchRProjects[[p]]@cellColData$timepoint_clone[ArchRProjects[[p]]@cellColData$timepoint=="REL"], replace = F)
    
    
    #ArchRProjects[[p]]<-addGroupCoverages(ArchRProjects[[p]], "timepoint_clone")
    
    peak_mtx_collapsed<-as.data.frame(matrix(nrow=nrow(peak_mtx), ncol=0))
    
    for(c in unique(ArchRProjects[[p]]@cellColData$timepoint_clone)){
      if(c %in% c("Unk_DX","Unk_REL")){next}
      
      cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clone==c]
      
      if(length(cells_of_interest)<50){next}
      
      cells_of_interest<-sample(cells_of_interest, 200, replace = T)
      
      pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
      colnames(pseudobulk)<-c(c)
      
      peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
    }
    
    if(p=="SU484"){
      c<-"02_DX"
      cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clone==c]
      pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
      colnames(pseudobulk)<-c(c)
      peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
    }
    
    peak_mtx_collapsed<-edgeR::cpm(peak_mtx_collapsed)
    
    sample_corr<-cor(peak_mtx_collapsed)
    
    cor_set[[i]]<-sample_corr
  }
  
  
  # generate actual correlations
  ArchRProjects[[p]]@cellColData$timepoint_clone<-paste(ArchRProjects[[p]]@cellColData$mito_clone,
                                                        ArchRProjects[[p]]@cellColData$timepoint, sep="_")
  
  #ArchRProjects[[p]]<-addGroupCoverages(ArchRProjects[[p]], "timepoint_clone")
  
  peak_mtx_collapsed<-as.data.frame(matrix(nrow=nrow(peak_mtx), ncol=0))
  
  for(c in unique(ArchRProjects[[p]]@cellColData$timepoint_clone)){
    if(c %in% c("Unk_DX","Unk_REL")){next}
    
    cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clone==c]
    
    if(length(cells_of_interest)<50){next}
    
    cells_of_interest<-sample(cells_of_interest, 200, replace = T)
    
    pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
    colnames(pseudobulk)<-c(c)
    
    peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
  }
  
  if(p=="SU484"){
    c<-"02_DX"
    cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clone==c]
    pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
    colnames(pseudobulk)<-c(c)
    peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
  }
  
  peak_mtx_collapsed<-edgeR::cpm(peak_mtx_collapsed)
  
  sample_corr<-cor(peak_mtx_collapsed)
  
  sample_corr_sub<-sample_corr[intersect(rownames(sample_corr),rownames(cor_set[[1]])),
                           intersect(colnames(sample_corr),colnames(cor_set[[1]]))]
  
  significance_mtx<-as.data.frame(matrix(0.5,nrow=nrow(sample_corr),ncol=ncol(sample_corr)))
  rownames(significance_mtx)<-rownames(sample_corr)
  colnames(significance_mtx)<-colnames(sample_corr)
  
  DX_cor<-sample_corr_sub[grep("DX",rownames(sample_corr_sub)),grep("DX",colnames(sample_corr_sub))]
  
  for(i in rownames(DX_cor)){
    for(j in colnames(DX_cor)){
      if(i!=j){
        message(paste(i,j))
        distribution<-sapply(1:100, function(x) cor_set[[x]][i,j])
        percentile <- ecdf(distribution)
        message(percentile(sample_corr_sub[i,j]))
        significance_mtx[i,j]<-percentile(sample_corr_sub[i,j])
      }
    }
  }
  
  REL_cor<-sample_corr_sub[grep("REL",rownames(sample_corr_sub)),grep("REL",colnames(sample_corr_sub))]
  
  for(i in rownames(REL_cor)){
    for(j in colnames(REL_cor)){
      if(i!=j){
        message(paste(i,j))
        distribution<-sapply(1:100, function(x) cor_set[[x]][i,j])
        percentile <- ecdf(distribution)
        message(percentile(sample_corr_sub[i,j]))
        significance_mtx[i,j]<-percentile(sample_corr_sub[i,j])
      }
    }
  }
  
  
  
  
  annotations<-data.frame(clone=strsplit(colnames(sample_corr),"_") %>% sapply('[',1),
                          timepoint=strsplit(colnames(sample_corr),"_") %>% sapply('[',2))
  clone_colors<-ArchRPalettes$stallion[1:length(unique(annotations$clone))]
  names(clone_colors)<-sort(unique(annotations$clone))
  timepoint_colors<-c("darkblue","darkred")
  names(timepoint_colors)<-c("DX","REL")
  anno_colors<-list(clone=clone_colors,
                    timepoint=timepoint_colors)
  
  ht<-Heatmap(sample_corr, 
              col=viridis(10),
              #col=colorRamp2(seq(0.5,1,length.out = 10), viridis(10)),
              show_row_names=FALSE, 
              show_column_names=T, 
              #show_column_dend = FALSE,
              top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"), 
              #left_annotation = rowAnnotation(df=annotations_mod, col=anno_colors, show_annotation_name=F),
              heatmap_legend_param = list(title = "Pearson"),
              column_title=paste(p,"Subclone Epigenetic Similarity"),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(significance_mtx[i, j] < 0.05) {
                  grid.text("*", x, y)
                }})
  
  
  pdf(paste0("outputs/scATAC_plots/mito_subclone_accessibility/",p,"_mito_subclone_similarity_significance_labeled.pdf"), width = 5, height = 4.5, useDingbats = FALSE)
  print(ht)
  dev.off()
}





# epigenetic difference between clones with differing degrees of difference with founder clone
degrees_difference<-list(SU142 = c(0,1,1),
                         SU360 = c(0,1,1,2),
                         SU484 = c(0,1,1,2,1,2,1),
                         SU892 = c(2,0,1,1))
origin_clones<-list(SU142 = c("00_DX"),
                    SU360 = c("00_DX"),
                    SU484 = c("00_DX"),
                    SU892 = c("01_DX"))

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
  
  
  
  #ArchRProjects[[p]]<-addGroupCoverages(ArchRProjects[[p]], "timepoint_clone")
  
  
  peak_mtx<-getMatrixFromProject(ArchRProj = ArchRProjects[[p]], useMatrix = "PeakMatrix")
  
  peak_mtx_collapsed<-as.data.frame(matrix(nrow=nrow(peak_mtx), ncol=0))
  
  for(c in unique(ArchRProjects[[p]]@cellColData$timepoint_clone)){
    if(c %in% c("Unk_DX","Unk_REL")){next}
    
    cells_of_interest<-rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint_clone==c]
    
    if(length(cells_of_interest)<100){next}
    
    pseudobulk<-data.frame(x=rowSums(peak_mtx[,cells_of_interest]@assays@data$PeakMatrix))
    colnames(pseudobulk)<-c(c)
    
    peak_mtx_collapsed<-cbind(peak_mtx_collapsed, pseudobulk)
  }
  
  peak_mtx_collapsed<-edgeR::cpm(peak_mtx_collapsed)
  
  sample_corr<-cor(peak_mtx_collapsed)
  
  degrees<-degrees_difference[[p]]
  
  origin_clone<-origin_clones[[p]]
  
  plt_input<-data.frame(similarity=c(1),
                        degrees_difference=c(0))
  for(clone in 1:length(degrees)){
    curr_clone<-paste0("0",as.character(clone-1),"_DX")
    if(curr_clone==origin_clone){next}
    if(!(curr_clone%in%colnames(sample_corr))){next}
    
    plt_input<-rbind(plt_input,
                     data.frame(similarity=c(sample_corr[origin_clone,curr_clone]),
                                degrees_difference=c(degrees[clone])))
  }
  
  plt_input$degrees_difference<-as.factor(plt_input$degrees_difference)
  
  plt<-ggplot(plt_input, aes(x=degrees_difference, y=similarity)) +
    geom_jitter(width = 0.2) +
    theme_classic() +
    ggtitle(p)
  print(plt)
}




# SU360 relapse accessibility signature comparison between clones 

p<-"SU360"

clone1_rel <- getMarkerFeatures(
  ArchRProj = ArchRProjects[[p]], 
  useMatrix = "GeneScoreMatrix",
  groupBy = "timepoint_clone",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "01_REL",
  bgdGroups = "01_DX"
)

clone1_rel <- getMarkers(clone1_rel, cutOff = "abs(Log2FC) >= 0")[["01_REL"]]
clone1_rel<-as.data.frame(clone1_rel) 
#clone1_rel$region<-paste(clone1_rel$seqnames, clone1_rel$start, clone1_rel$end, sep="_")

clone2_rel <- getMarkerFeatures(
  ArchRProj = ArchRProjects[[p]], 
  useMatrix = "GeneScoreMatrix",
  groupBy = "timepoint_clone",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "02_REL",
  bgdGroups = "02_DX"
)

clone2_rel <- getMarkers(clone2_rel, cutOff = "abs(Log2FC) >= 0")[["02_REL"]]
clone2_rel<-as.data.frame(clone2_rel) 
#clone2_rel$region<-paste(clone2_rel$seqnames, clone2_rel$start, clone2_rel$end, sep="_")


plot_input<-data.frame(region=clone1_rel$name,
                       clone_1_lfc=clone1_rel$Log2FC,
                       clone_2_lfc=clone2_rel$Log2FC[match(clone1_rel$name, clone2_rel$name)])

plt<-ggplot(plot_input, aes(x=clone_1_lfc, y=clone_2_lfc)) +
  ggrastr::rasterise(geom_point(size=0.5), dpi=200) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5)) +
  scale_x_continuous(limits = c(-2,2)) +
  scale_y_continuous(limits = c(-2,2)) +
  stat_cor(method = "pearson", label.sep = "\n", size=4) +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  xlab("Clone 1 Relapse vs Diagnosis Log Fold Change") +
  ylab("Clone 2 Relapse vs Diagnosis Log Fold Change") +
  ggtitle("SU360 Clone 1 vs Clone 2\nRelapse Fold Change")

pdf(paste0("outputs/scATAC_plots/mito_subclone_accessibility/SU360_clone1_vs_clone2_REL_vs_DX_accessibility_dotplot.pdf"), width = 4, height = 4, useDingbats = FALSE)
print(plt)
dev.off()

# SU484 differential accessibility between clones that survive and die at relapse 

p<-"SU484"


ArchRProjects[[p]]@cellColData$survival_at_relapse<-rep(NA, nrow(ArchRProjects[[p]]@cellColData))
ArchRProjects[[p]]@cellColData$survival_at_relapse[ArchRProjects[[p]]@cellColData$timepoint_clone%in%c("01_DX", "03_DX","04_DX", "05_DX","06_DX")]<-"eliminated"
ArchRProjects[[p]]@cellColData$survival_at_relapse[ArchRProjects[[p]]@cellColData$timepoint_clone%in%c("00_DX")]<-"survived"


diff_res <- getMarkerFeatures(
  ArchRProj = ArchRProjects[[p]], 
  #useGroups = c("00_DX","01_DX","02_DX", "03_DX","04_DX", "05_DX","06_DX"),
  useGroups = c("00_REL","02_REL"),
  useMatrix = "GeneScoreMatrix",
  groupBy = "timepoint_clone",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)")
)

markerList <- getMarkers(diff_res, cutOff = "FDR <= 1 & Log2FC >= 0")
markerList<-as.data.frame(markerList)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = diff_res, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  plotLog2FC = T, nLabel = 30,
  transpose = F
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

markerPlot(seMarker = diff_res, name = "00_DX", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "Volcano")


# tracks of select regions across clones

gene  <- c("ERG")
HOXA_locus<-GRanges(seqnames = c("chr7"), ranges = IRanges(start=c(27080000), end=c(27220000)))
HOXB_locus<-GRanges(seqnames = c("chr17"), ranges = IRanges(start=c(48520000), end=c(48640000)))
HIST_locus<-GRanges(seqnames = c("chr6"), ranges = IRanges(start=c(26000000), end=c(26300000)))
CDK6_locus<-GRanges(seqnames = c("chr7"), ranges = IRanges(start=c(92489398), end=c(93088969)))

plt <- plotBrowserTrack(
  ArchRProj = ArchRProjects[[p]], 
  groupBy = "timepoint_clone",
  #region = HIST_locus,
  geneSymbol = gene, minCells = 0,
  useGroups = c("00_DX","00_REL","01_DX","02_DX","02_REL","03_DX","04_DX","05_DX"),
  upstream = 100000,
  downstream = 100000, loops = NULL
)

grid::grid.newpage()
grid::grid.draw(plt)
grid::grid.draw(plt[[gene]])




ArchRProjects[[p]] <- addCoAccessibility(
  ArchRProj = ArchRProjects[[p]], 
  cellsToUse = rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint=="DX"],
  reducedDims = "IterativeLSI"
)

dx_coaccessibility<-getCoAccessibility(ArchRProjects[[p]])

ArchRProjects[[p]] <- addCoAccessibility(
  ArchRProj = ArchRProjects[[p]], 
  cellsToUse = rownames(ArchRProjects[[p]]@cellColData)[ArchRProjects[[p]]@cellColData$timepoint=="REL"],
  reducedDims = "IterativeLSI"
)

rel_coaccessibility<-getCoAccessibility(ArchRProjects[[p]])


gene<-"MEIS1"

plt1 <- plotBrowserTrack(
  ArchRProj = ArchRProjects[[p]], 
  groupBy = "Sample",
  #region = CDK6_locus,
  geneSymbol = gene,
  #useGroups = c("00_DX","00_REL","01_DX","01_REL","02_DX","02_REL","03_DX","03_REL","04_DX","04_REL","05_DX","05_REL"),
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  loops = dx_coaccessibility,
  upstream = 500000,
  downstream = 500000
)

plt2 <- plotBrowserTrack(
  ArchRProj = ArchRProjects[[p]], 
  groupBy = "Sample",
  #region = CDK6_locus,
  geneSymbol = gene,
  #useGroups = c("00_DX","00_REL","01_DX","01_REL","02_DX","02_REL","03_DX","03_REL","04_DX","04_REL","05_DX","05_REL"),
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  loops = rel_coaccessibility,
  upstream = 500000,
  downstream = 500000, 
)

grid::grid.newpage()
grid::grid.draw(plt1[[gene]])

grid::grid.newpage()
grid::grid.draw(plt2[[gene]])

#ggarrange(plt1[[gene]],plt2[[gene]], ncol=1)



length(intersect(paste(dx_coaccessibility$CoAccessibility@seqnames,dx_coaccessibility$CoAccessibility@ranges@start,dx_coaccessibility$CoAccessibility@ranges@width,sep="_"),
                 paste(rel_coaccessibility$CoAccessibility@seqnames,rel_coaccessibility$CoAccessibility@ranges@start,rel_coaccessibility$CoAccessibility@ranges@width,sep="_")))
