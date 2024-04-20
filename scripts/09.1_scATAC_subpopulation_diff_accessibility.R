# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Differential accessibility between DX and REL within different AML subpopulations
# subpopulations defined by projection and nearest neighbor classification
#


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

# do differential accessiblity for all cells between timepoints
for(p in unique(scATAC_ref$patient)){
  
  message(p)
  
  ArchRProjects[[p]]@cellColData$timepoint<-rep("DX", nrow(ArchRProjects[[p]]@cellColData))
  ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
  
  diff_res <- getMarkerFeatures(
    ArchRProj = ArchRProjects[[p]], 
    useMatrix = "GeneScoreMatrix",
    groupBy = "timepoint",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "REL",
    bgdGroups = "DX"
  )
    
    markerList <- getMarkers(diff_res, cutOff = "abs(Log2FC) >= 0")[["REL"]]
    markerList<-as.data.frame(markerList) 
    
    saveRDS(markerList, paste0("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes/",p,"_bulk_REL_vs_DX_differential_gene_accessibility.rds"))
}

# do differential accessibility for each subpopulation
for(p in unique(scATAC_ref$patient)){
  
  message(p)
  
  # label classification celltypes
  classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  ArchRProjects[[p]]@cellColData$projection_classification<-classifications$AML_classification[match(rownames(ArchRProjects[[p]]@cellColData),rownames(classifications))]
  ArchRProjects[[p]]@cellColData$projection_classification[is.na(ArchRProjects[[p]]@cellColData$projection_classification)]<-"Unk"
  
  ArchRProjects[[p]]@cellColData$timepoint<-rep("DX", nrow(ArchRProjects[[p]]@cellColData))
  ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
  
  ArchRProjects[[p]]@cellColData$timepoint_classification<-paste(ArchRProjects[[p]]@cellColData$projection_classification,
                                                                 ArchRProjects[[p]]@cellColData$timepoint, sep="_")
  
  
  celltypes<-unique(ArchRProjects[[p]]@cellColData$projection_classification)
  
  for(c in celltypes){
    
    if(c=="Unk"){next}
    
    message(c)
    
    if(sum(ArchRProjects[[p]]@cellColData$projection_classification==c)<100){
      message("less than 100 cells, skipping")
      next
    }
    
    diff_res <- getMarkerFeatures(
      ArchRProj = ArchRProjects[[p]], 
      useMatrix = "GeneScoreMatrix",
      groupBy = "timepoint_classification",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = paste0(c,"_REL"),
      bgdGroups = paste0(c,"_DX")
    )
    
    markerList <- getMarkers(diff_res, cutOff = "abs(Log2FC) >= 0")[[paste0(c,"_REL")]]
    markerList<-as.data.frame(markerList) 
    
    c<-gsub("/","-",c)
    
    saveRDS(markerList, paste0("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes/",p,"_",c,"_REL_vs_DX_differential_gene_accessibility.rds"))
  }
}


# plot heatmap of the differential genes for each subpopulation

files<-list.files("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes")
lfc_df<-as.data.frame(matrix(nrow=nrow(readRDS("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes/SU142_CMP-LMPP-like_REL_vs_DX_differential_gene_accessibility.rds")),
                             ncol=0))
rownames(lfc_df)<-readRDS("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes/SU142_CMP-LMPP-like_REL_vs_DX_differential_gene_accessibility.rds")$name

for(f in files){
  p<-strsplit(f, "_") %>% sapply('[',1)
  c<-strsplit(f, "_") %>% sapply('[',2)
  
  temp_df<-readRDS(paste0("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes/",f))
  temp_df_1<-temp_df[,c("Log2FC"),drop=F]
  rownames(temp_df_1)<-temp_df$name
  temp_df_1<-temp_df_1[rownames(lfc_df),,drop=F]
  colnames(temp_df_1)<-c(paste(p,c,sep="_"))
  
  lfc_df<-cbind(lfc_df,
                temp_df_1)
}
lfc_df<-lfc_df[complete.cases(lfc_df),]
lfc_df<-lfc_df[apply(lfc_df,1, sd)>0,]


celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP-LMPP-like","Mono-DC-Baso-like","GMP-like","CLP-like","B-Plasma-like","Erythroid-like","T-NK-like"))

# similarity analysis of differential chromatin accessibility

pca_res<-as.data.frame(prcomp(lfc_df)$rotation)
pca_res$patient<-strsplit(rownames(pca_res), "_") %>% sapply('[',1)
pca_res$celltype<-strsplit(rownames(pca_res), "_") %>% sapply('[',2)

ggplot(pca_res, aes(x=PC1, y=PC3, color=celltype, shape=patient)) +
  geom_point() +
  theme_classic()


sample_corr<-cor(lfc_df)

annotations<-data.frame(patient=strsplit(colnames(sample_corr),"_") %>% sapply('[',1),
                        celltype=strsplit(colnames(sample_corr),"_") %>% sapply('[',2))
celltype_colors<-c(celltype_pal,"grey")
names(celltype_colors)<-c(names(celltype_colors)[1:(length(celltype_colors)-1)],"bulk")
patient_colors<-get_palette("ucscgb",length(unique(annotations$patient)))
names(patient_colors)<-unique(annotations$patient)
anno_colors<-list(celltype=celltype_colors,
                  patient=patient_colors)

ht<-Heatmap(sample_corr, 
            #col=viridis(10),
            col=colorRamp2(seq(0,0.8,length.out = 10), viridis(10)),
            show_row_names=FALSE, 
            show_column_names=T, 
            #show_column_dend = FALSE,
            top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"), 
            #left_annotation = rowAnnotation(df=annotations_mod, col=anno_colors, show_annotation_name=F),
            heatmap_legend_param = list(title = "Pearson"),
            column_title="Subpopulation-specific DX vs REL LFC Similarity",
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8))


print(ht)

pdf("outputs/scATAC_plots/scATAC_subpoopulation_differential_accessibility/all_subpopulation_fold_change_global_similarity.pdf", width = 6.5, height = 6, useDingbats = FALSE)
print(ht)
dev.off()






# 
# fdr_df<-as.data.frame(matrix(nrow=nrow(readRDS("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes/SU142_CMP-LMPP-like_REL_vs_DX_differential_gene_accessibility.rds")),
#                              ncol=0))
# rownames(fdr_df)<-readRDS("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes/SU142_CMP-LMPP-like_REL_vs_DX_differential_gene_accessibility.rds")$name
# 
# for(f in files){
#   p<-strsplit(f, "_") %>% sapply('[',1)
#   c<-strsplit(f, "_") %>% sapply('[',2)
#   
#   temp_df<-readRDS(paste0("outputs/scATAC/AML_scATAC_diff_accessibility_closest_celltypes/",f))
#   temp_df_1<-temp_df[,c("FDR"),drop=F]
#   rownames(temp_df_1)<-temp_df$name
#   temp_df_1<-temp_df_1[rownames(fdr_df),,drop=F]
#   colnames(temp_df_1)<-c(paste(p,c,sep="_"))
#   
#   fdr_df<-cbind(fdr_df,
#                 temp_df_1)
# }
# fdr_df<-fdr_df[complete.cases(fdr_df),]
# 
# significant_features_by_fdr<-rownames(fdr_df[apply(fdr_df,1, min)<0.01,])

annotations<-data.frame(patient=strsplit(files, "_") %>% sapply('[',1),
                        celltype=strsplit(files, "_") %>% sapply('[',2))
annotations$name<-paste(annotations$patient,annotations$celltype,sep="_")

# gene_disp<-get_variable_gene(t(lfc_df))
# 
# heatmap_input<-lfc_df[rownames(gene_disp[order(gene_disp$dispersion_norm,decreasing = T),])[1:1000],]
# 
# group = kmeans(heatmap_input, centers = 5)$cluster
# 
# Heatmap(heatmap_input, 
#         cluster_rows = cluster_within_group(t(heatmap_input), group),
#         #cluster_rows = T,
#         col=colorRamp2(seq(-1,1,length.out = 10), viridis(10)),
#         show_row_names=FALSE, 
#         show_column_names=T,
#         #heatmap_legend_param = list(title = "Row Z-Scores"),
#         #column_title="Stable Samples Differential Peak Accessibility",
#         #top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"),
#         row_names_gp = gpar(fontsize = 8),
#         column_names_gp = gpar(fontsize = 8))



## ID featires that are patient specific

pt_order<-as.factor(unique(annotations$patient))
assignments<-annotations[,c("patient"),drop=FALSE]
rownames(assignments)<-annotations$name

feature_assignments<-get_feature_assignments(lfc_df, assignments, pt_order)
feature_assignments<-feature_assignments[complete.cases(feature_assignments),]

significant_features<-c()
plt_genes<-c()
max_genes<-500

for(c in unique(assignments$patient)){
  message(c)

  grp1<-rownames(assignments)[assignments$patient==c]
  grp2<-rownames(assignments)[assignments$patient!=c]

  lfc<-apply(lfc_df[feature_assignments$feature[feature_assignments$assignment==c],],1,function(x) mean(x[grp1])-mean(x[grp2]))

  feature_info<-data.frame(feature=rownames(lfc_df[feature_assignments$feature[feature_assignments$assignment==c],]),
                           lfc=lfc)


  if(nrow(feature_info)>max_genes){feature_info<-feature_info[rev(order(feature_info$lfc)),][1:max_genes,]}

  message(length(feature_info$feature))

  significant_features<-c(significant_features, feature_info$feature)
  plt_genes<-c(plt_genes,feature_info$feature[order(feature_info$lfc, decreasing = T)][1:10])
}

plot_input<-lfc_df[significant_features,rownames(assignments)]
#plot_input<-as.data.frame(t(scale(t(plot_input))))
plot_input<-reorder_dataframe(plot_input, assignments, pt_order)


feature_groups<-feature_assignments$assignment
names(feature_groups)<-feature_assignments$feature
feature_groups<-feature_groups[rownames(plot_input)]

plt_genes_idx<-match(plt_genes, rownames(plot_input))


if(T){
  plt_annotations<-column_to_rownames(annotations, "name")
  plt_annotations<-plt_annotations[colnames(plot_input),,drop=FALSE]
  plt_annotations$patient<-as.factor(plt_annotations$patient)
  patient_colors<-get_palette("ucscgb", length(unique(plt_annotations$patient)))
  names(patient_colors)<-unique(plt_annotations$patient)
  celltype_colors<-c(celltype_pal,"grey")
  names(celltype_colors)<-c(names(celltype_colors)[1:(length(celltype_colors)-1)],"bulk")
  anno_colors<-list(celltype=celltype_colors,
                    patient=patient_colors)
  
  col_anno<-data.frame(patient=feature_groups)
  rownames(col_anno)<-names(feature_groups)
  col_anno<-col_anno[rownames(plot_input),,drop=FALSE]
  col_anno_colors<-list(patient=patient_colors)
  
}

Heatmap(plot_input, 
        #cluster_rows = cluster_within_group(t(plot_input), feature_groups),
        cluster_rows = F,
        cluster_columns = F,
        col=colorRamp2(seq(-0.5,0.5,length.out = 10), viridis(10)),
        show_row_names=F, 
        show_column_names=T,
        heatmap_legend_param = list(title = "log2(REL-DX Fold Change)"),
        #column_title="Stable Samples Differential Peak Accessibility",
        top_annotation=HeatmapAnnotation(df=plt_annotations, col=anno_colors, annotation_name_side = "left"),
        right_annotation = rowAnnotation(genes = anno_mark(at = plt_genes_idx, 
                                                           labels = plt_genes, 
                                                           side="right",
                                                           labels_gp = gpar(fontsize = 8),
                                                           link_width =unit(5, "mm"),
                                                           padding = 0.7)),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))





## ID featires that are celltype specific

lfc_df_sub<-lfc_df[,grep("bulk",colnames(plot_input),value=T, invert = T)]

assignments<-annotations[,c("celltype"),drop=FALSE]
rownames(assignments)<-annotations$name
assignments<-assignments[assignments$celltype%in%c("CMP-LMPP-like","GMP-like","Mono-DC-Baso-like"),,drop=F]

assignments<-assignments[assignments$celltype!="bulk",,drop=F]
celltype_order<-as.factor(unique(assignments$celltype))

feature_assignments<-get_feature_assignments(lfc_df, assignments, celltype_order)
feature_assignments<-feature_assignments[complete.cases(feature_assignments),]

significant_features<-c()
plt_genes<-c()
max_genes<-100

for(c in unique(assignments$celltype)){
  message(c)
  
  grp1<-rownames(assignments)[assignments$celltype==c]
  grp2<-rownames(assignments)[assignments$celltype!=c]
  
  lfc<-apply(lfc_df_sub[feature_assignments$feature[feature_assignments$assignment==c],],1,function(x) mean(as.numeric(x[grp1]))-mean(as.numeric(x[grp2])))
  
  feature_info<-data.frame(feature=rownames(lfc_df_sub[feature_assignments$feature[feature_assignments$assignment==c],]),
                           lfc=lfc)
  
  
  if(nrow(feature_info)>max_genes){feature_info<-feature_info[rev(order(feature_info$lfc)),][1:max_genes,]}
  
  message(length(feature_info$feature))
  
  significant_features<-c(significant_features, feature_info$feature)
  plt_genes<-c(plt_genes,feature_info$feature[order(feature_info$lfc, decreasing = T)][1:10])
}

plot_input<-lfc_df_sub[significant_features,rownames(assignments)]
#plot_input<-as.data.frame(t(scale(t(plot_input))))
plot_input<-reorder_dataframe(plot_input, assignments, celltype_order)


feature_groups<-feature_assignments$assignment
names(feature_groups)<-feature_assignments$feature
feature_groups<-feature_groups[rownames(plot_input)]

plt_genes_idx<-match(plt_genes, rownames(plot_input))

if(T){
  plt_annotations<-column_to_rownames(annotations, "name")
  plt_annotations<-plt_annotations[colnames(plot_input),,drop=FALSE]
  plt_annotations$patient<-as.factor(plt_annotations$patient)
  patient_colors<-get_palette("ucscgb", length(unique(plt_annotations$patient)))
  names(patient_colors)<-unique(plt_annotations$patient)
  celltype_colors<-c(celltype_pal,"grey")
  names(celltype_colors)<-c(names(celltype_colors)[1:(length(celltype_colors)-1)],"bulk")
  anno_colors<-list(celltype=celltype_colors,
                    patient=patient_colors)
  
  col_anno<-data.frame(patient=feature_groups)
  rownames(col_anno)<-names(feature_groups)
  col_anno<-col_anno[rownames(plot_input),,drop=FALSE]
  col_anno_colors<-list(patient=patient_colors)
  
}

Heatmap(plot_input, 
        #cluster_rows = cluster_within_group(t(plot_input), feature_groups),
        cluster_rows = F,
        cluster_columns = F,
        col=colorRamp2(seq(-0.5,0.5,length.out = 10), viridis(10)),
        show_row_names=F, 
        show_column_names=T,
        heatmap_legend_param = list(title = "log2(REL-DX Fold Change)"),
        #column_title="Stable Samples Differential Peak Accessibility",
        top_annotation=HeatmapAnnotation(df=plt_annotations, col=anno_colors, annotation_name_side = "left"),
        right_annotation = rowAnnotation(genes = anno_mark(at = plt_genes_idx, 
                                                           labels = plt_genes, 
                                                           side="right",
                                                           labels_gp = gpar(fontsize = 8),
                                                           link_width =unit(5, "mm"),
                                                           padding = 0.7)),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))















#deseq_Results<-readRDS("outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")
#top_genes<-rownames(relapse_signature)[relapse_signature$padj<0.05&abs(relapse_signature$log2FoldChange)>0.5] %>% unique()

# top_regions<-deseq_Results[deseq_Results$log2FoldChange>0&deseq_Results$padj<0.05,]
# top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:100]
# bottom_regions<-deseq_Results[deseq_Results$log2FoldChange<0&deseq_Results$padj<0.05,]
# bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:100]

# plot_input<-lfc_df[significant_features_by_fdr,]
# plot_input<-lfc_df[c(intersect(c(top_regions), rownames(lfc_df)),
#                      intersect(bottom_regions, rownames(lfc_df))),]
plot_input<-lfc_df[sample(rownames(lfc_df), 2000), grep("SU360",colnames(lfc_df), value = T)]
plot_input<-plot_input[order(rowSums(plot_input), decreasing = T),]
plot_input<-lfc_df[sample(rownames(lfc_df), 2000),]

#top_regions<-rownames(lfc_df)[order(lfc_df$SU360_bulk, decreasing = T)][1:200]
#top_regions<-rownames(lfc_df)[order(rowSums(lfc_df[,grep("SU360",colnames(lfc_df))]), decreasing = T)][1:200]
#plot_input<-lfc_df[top_regions,grep("SU360",colnames(lfc_df), value = T)]

plot_input<-plot_input[,grep("GMP|LMPP",colnames(plot_input),value=T, invert = F)]

feature_groups<-feature_assignments$assignment
names(feature_groups)<-feature_assignments$feature
feature_groups<-feature_groups[rownames(plot_input)]

if(T){
  plt_annotations<-column_to_rownames(annotations, "name")
  plt_annotations<-plt_annotations[colnames(plot_input),,drop=FALSE]
  plt_annotations$patient<-as.factor(plt_annotations$patient)
  patient_colors<-get_palette("ucscgb", length(unique(plt_annotations$patient)))
  names(patient_colors)<-unique(plt_annotations$patient)
  celltype_colors<-ArchR::paletteDiscrete(c("HSC-like","CMP-LMPP-like","Mono-DC-Baso-like","GMP-like","CLP-like","B-Plasma-like","Erythroid-like","T-NK-like","bulk"))
  anno_colors<-list(patient=patient_colors, celltype=celltype_colors)
  
  col_anno<-data.frame(patient=feature_groups)
  rownames(col_anno)<-names(feature_groups)
  col_anno<-col_anno[rownames(plot_input),,drop=FALSE]
  col_anno_colors<-list(patient=patient_colors)
  
}

Heatmap(plot_input, 
        #cluster_rows = cluster_within_group(t(plot_input), feature_groups),
        cluster_rows = T,
        cluster_columns = F,
        col=colorRamp2(seq(-0.5,0.5,length.out = 10), viridis(10)),
        show_row_names=F, 
        show_column_names=T,
        heatmap_legend_param = list(title = "log2(REL-DX Fold Change)"),
        #column_title="Stable Samples Differential Peak Accessibility",
        top_annotation=HeatmapAnnotation(df=plt_annotations, col=anno_colors, annotation_name_side = "left"),
        #right_annotation=rowAnnotation(df=col_anno, col=col_anno_colors),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))






# motif analysis

for(p in unique(scATAC_ref$patient)){
  
  message(p)
  
  # label classification celltypes
  classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  ArchRProjects[[p]]@cellColData$projection_classification<-classifications$AML_classification[match(rownames(ArchRProjects[[p]]@cellColData),rownames(classifications))]
  ArchRProjects[[p]]@cellColData$projection_classification[is.na(ArchRProjects[[p]]@cellColData$projection_classification)]<-"Unk"
  
  ArchRProjects[[p]]@cellColData$timepoint<-rep("DX", nrow(ArchRProjects[[p]]@cellColData))
  ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
  
  ArchRProjects[[p]]@cellColData$timepoint_classification<-paste(ArchRProjects[[p]]@cellColData$projection_classification,
                                                                 ArchRProjects[[p]]@cellColData$timepoint, sep="_")
  
  
  celltypes<-unique(ArchRProjects[[p]]@cellColData$projection_classification)
  
  for(c in celltypes){
    
    if(c=="Unk"){next}
    
    message(c)
    
    if(sum(ArchRProjects[[p]]@cellColData$projection_classification==c)<100){
      message("less than 100 cells, skipping")
      next
    }
    
    diff_res <- getMarkerFeatures(
      ArchRProj = ArchRProjects[[p]], 
      useMatrix = "PeakMatrix",
      groupBy = "timepoint_classification",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = paste0(c,"_REL"),
      bgdGroups = paste0(c,"_DX")
    )
    
    motifsUp <- peakAnnoEnrichment(
      seMarker = diff_res,
      ArchRProj = ArchRProjects[[p]],
      peakAnnotation = "Motif",
      cutOff = "FDR <= 0.1 & Log2FC >0"
    )
    
    motifsDn <- peakAnnoEnrichment(
      seMarker = diff_res,
      ArchRProj = ArchRProjects[[p]],
      peakAnnotation = "Motif",
      cutOff = "FDR <= 0.1 & Log2FC < 0"
    )
    
    motifs_res <- data.frame(TF = rownames(motifsUp), 
                             mlog10Padj_up = assay(motifsUp)[,1],
                             mlog10Padj_dn = assay(motifsDn)[,1])
    
    c<-gsub("/","-",c)
    
    saveRDS(motifs_res, paste0("outputs/scATAC/AML_scATAC_diff_motifs_closest_celltypes/",p,"_",c,"_REL_vs_DX_differential_motif_accessibility.rds"))
  }
}


# generate motif matrix
files<-list.files("outputs/scATAC/AML_scATAC_diff_motifs_closest_celltypes")

motif_df_up<-as.data.frame(matrix(nrow=nrow(readRDS("outputs/scATAC/AML_scATAC_diff_motifs_closest_celltypes/SU142_CMP-LMPP-like_REL_vs_DX_differential_motif_accessibility.rds")), ncol=0))
rownames(motif_df_up)<-readRDS("outputs/scATAC/AML_scATAC_diff_motifs_closest_celltypes/SU142_CMP-LMPP-like_REL_vs_DX_differential_motif_accessibility.rds")$TF
for(f in files){
  p<-strsplit(f, "_") %>% sapply('[',1)
  c<-strsplit(f, "_") %>% sapply('[',2)
  
  temp_df<-readRDS(paste0("outputs/scATAC/AML_scATAC_diff_motifs_closest_celltypes/",f))
  temp_df_1<-temp_df[,c("mlog10Padj_up"),drop=F]
  rownames(temp_df_1)<-temp_df$TF
  temp_df_1<-temp_df_1[rownames(motif_df_up),,drop=F]
  colnames(temp_df_1)<-c(paste(p,c,sep="_"))
  
  motif_df_up<-cbind(motif_df_up,
                temp_df_1)
}

motif_df_dn<-as.data.frame(matrix(nrow=nrow(readRDS("outputs/scATAC/AML_scATAC_diff_motifs_closest_celltypes/SU142_CMP-LMPP-like_REL_vs_DX_differential_motif_accessibility.rds")), ncol=0))
rownames(motif_df_dn)<-readRDS("outputs/scATAC/AML_scATAC_diff_motifs_closest_celltypes/SU142_CMP-LMPP-like_REL_vs_DX_differential_motif_accessibility.rds")$TF
for(f in files){
  p<-strsplit(f, "_") %>% sapply('[',1)
  c<-strsplit(f, "_") %>% sapply('[',2)
  
  temp_df<-readRDS(paste0("outputs/scATAC/AML_scATAC_diff_motifs_closest_celltypes/",f))
  temp_df_1<-temp_df[,c("mlog10Padj_dn"),drop=F]
  rownames(temp_df_1)<-temp_df$TF
  temp_df_1<-temp_df_1[rownames(motif_df_dn),,drop=F]
  colnames(temp_df_1)<-c(paste(p,c,sep="_"))
  
  motif_df_dn<-cbind(motif_df_dn,
                  temp_df_1)
}

# motif_df<-motif_df[apply(motif_df,1, sd)>0,]
# 
# motif_df<-motif_df[apply(motif_df,1, max)>5,]
# 
# temp<-motif_df
# temp[temp<5]<-0
# temp[temp>5]<-1
# motif_df<-motif_df[rowSums(temp)>=2,]

tf_search<-c("GATA","CEBP","RUNX","HOX","TAL","BCL")

interesting_tfs<-grep(paste(tf_search, collapse = "|"), rownames(motif_df_up), value = T)

plot_input<-motif_df_up#-motif_df_dn

plot_input<-plot_input[,grep("SU142",colnames(plot_input))]

#plot_input<-plot_input[interesting_tfs,]

plot_input<-plot_input[apply(plot_input,1, sd)>0,]

plot_input<-plot_input[apply(plot_input,1, function(x) max(abs(x)))>2,]

plot_input<-plot_input[grep("FOX",rownames(plot_input),value=T,invert = T),]

annotations<-data.frame(patient=strsplit(files, "_") %>% sapply('[',1),
                        celltype=strsplit(files, "_") %>% sapply('[',2))
annotations$name<-paste(annotations$patient,annotations$celltype,sep="_")

if(T){
  plt_annotations<-column_to_rownames(annotations, "name")
  plt_annotations<-plt_annotations[colnames(plot_input),,drop=FALSE]
  plt_annotations$patient<-as.factor(plt_annotations$patient)
  patient_colors<-get_palette("ucscgb", length(unique(plt_annotations$patient)))
  names(patient_colors)<-unique(plt_annotations$patient)
  celltype_colors<-ArchR::paletteDiscrete(c("HSC-like","CMP-LMPP-like","Mono-DC-Baso-like","GMP-like","CLP-like","B-Plasma-like","Erythroid-like","T-NK-like"))
  anno_colors<-list(patient=patient_colors, celltype=celltype_colors)
}

#plot_input<-plot_input[,grep("SU142",colnames(plot_input),value=T)]

Heatmap(plot_input, 
        #cluster_rows = cluster_within_group(t(plot_input), feature_groups),
        cluster_rows = T,
        cluster_columns = T,
        col=colorRamp2(seq(-2,2,length.out = 10), viridis(10)),
        show_row_names=T, 
        show_column_names=T,
        #heatmap_legend_param = list(title = "Row Z-Scores"),
        #column_title="Stable Samples Differential Peak Accessibility",
        top_annotation=HeatmapAnnotation(df=plt_annotations, col=anno_colors, annotation_name_side = "left"),
        #right_annotation=rowAnnotation(df=col_anno, col=col_anno_colors),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8))
















for(p in unique(scATAC_ref$patient)){
  
  message(p)
  
  # label classification celltypes
  classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  ArchRProjects[[p]]@cellColData$projection_classification<-classifications$AML_classification[match(rownames(ArchRProjects[[p]]@cellColData),rownames(classifications))]
  ArchRProjects[[p]]@cellColData$projection_classification[is.na(ArchRProjects[[p]]@cellColData$projection_classification)]<-"Unk"
  
  
  ArchRProjects[[p]]@cellColData$timepoint<-rep("DX", nrow(ArchRProjects[[p]]@cellColData))
  ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
  
  ArchRProjects[[p]]@cellColData$timepoint_classification<-paste(ArchRProjects[[p]]@cellColData$projection_classification,
                                                                 ArchRProjects[[p]]@cellColData$timepoint, sep="_")
  
  cell_counts<-as.data.frame(table(ArchRProjects[[p]]@cellColData$timepoint_classification))
  
  diff_res <- getMarkerFeatures(
    ArchRProj = ArchRProjects[[p]], 
    useMatrix = "PeakMatrix",
    groupBy = "timepoint_classification",
    useGroups = as.character(cell_counts$Var1[cell_counts$Freq>100]),
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
  )
  
  
  enrichMotifs <- peakAnnoEnrichment(
    seMarker = diff_res,
    ArchRProj = ArchRProjects[[p]], 
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0"
  )
  
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  
}
