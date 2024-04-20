# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Plotting and analysis of scATAC mitochondrial data
#
# Some code borrowed from Caleb Lareau
#
# Workflow for this script is as follows:
# 1. 
# 2. 
# 3. 

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

library(ComplexHeatmap)
library(viridis)
library(Seurat)
library(stringr)
library(ggpubr)
library(ArchR)

files<-list.files("inputs/scATAC/mito_vars_data/") %>% grep("rds", ., value = T)

file_reference<-data.frame(patient=substr(files,0,5),
                           sample=sapply(strsplit(files, "-|_"),'[',1),
                           file=files)
file_reference$timepoint<-rep("DX",nrow(file_reference))
file_reference$timepoint[nchar(file_reference$sample)==6]<-"REL"


# manually curate mt mutations based on post-filtering plots
# manual_mutations<-list(
#   SU142=c("16215A>G", "310T>C", "933G>A", "13369T>C", "4429G>A",  "3496G>A", "1107A>G"),
#   SU332=c("13129C>T", "11719G>A", "8251G>A", "11928A>G", "3921C>A", "1719G>A", "1703C>T", "14581T>C", "4960C>T", "310T>C", "16131T>C", "9310T>C"),
#   #SU360=c("2471G>A", "6810G>A", "515A>G", "11838T>A", "2805A>T", "7775G>A", "6456G>A", "12871G>A", "1040T>C", "7849C>T", "2908T>C"),
#   SU360=c("2471G>A", "6810G>A", "11838T>A", "2805A>T", "7775G>A", "6456G>A", "12871G>A", "1040T>C", "7849C>T", "2908T>C"),
#   #SU484=c("12476G>A", "310T>C", "515A>G", "3244G>A", "7069T>C", "15713T>C", "3199T>C", "6472T>C"),
#   SU484=c("12476G>A", "310T>C", "3244G>A", "7069T>C", "15713T>C", "3199T>C", "6472T>C"),
#   SU892=c("310T>C", "15639T>G", "16189T>C"),
#   SU926=c("9210A>G", "12414T>C", "7061A>G", "310T>C", "11955A>G", "515A>G", "523A>C")
# )

manual_mutations<-list(
  SU142=c("16215A>G", "933G>A", "13369T>C", "4429G>A",  "3496G>A", "1107A>G"),
  SU332=c("13129C>T", "11719G>A", "8251G>A", "11928A>G", "3921C>A", "1719G>A", "1703C>T", "14581T>C", "4960C>T", "310T>C", "16131T>C", "9310T>C"),
  #SU360=c("2471G>A", "6810G>A", "515A>G", "11838T>A", "2805A>T", "7775G>A", "6456G>A", "12871G>A", "1040T>C", "7849C>T", "2908T>C"),
  SU360=c("2471G>A", "6810G>A", "11838T>A", "2805A>T", "7775G>A", "6456G>A", "12871G>A", "1040T>C", "7849C>T", "2908T>C"),
  #SU484=c("12476G>A", "515A>G", "3244G>A", "7069T>C", "15713T>C", "3199T>C", "6472T>C"),
  SU484=c("12476G>A", "3244G>A", "7069T>C", "15713T>C", "3199T>C", "6472T>C"),
  SU892=c("15639T>G", "16189T>C"),
  SU926=c("9210A>G", "12414T>C", "7061A>G", "11955A>G", "515A>G", "523A>C")
)

manual_resolution<-list(SU142=0.25,
                        SU332=0.25,
                        SU360=0.05,
                        SU484=0.1,
                        SU892=0.01,
                        SU926=0.1)


for(p in unique(file_reference$patient)){
  message(p)
  
  dx_mito_vars<-readRDS(paste0("inputs/scATAC/mito_vars_data/",file_reference$file[file_reference$patient==p&file_reference$timepoint=="DX"]))
  rel_mito_vars<-readRDS(paste0("inputs/scATAC/mito_vars_data/",file_reference$file[file_reference$patient==p&file_reference$timepoint=="REL"]))
  
  colnames(dx_mito_vars)<-paste0("DX_",colnames(dx_mito_vars))
  colnames(rel_mito_vars)<-paste0("REL_",colnames(rel_mito_vars))
  
  mut_df<-data.frame(n_cells_conf_detected=data.frame(rowData(dx_mito_vars))$n_samples_conf_detected+data.frame(rowData(rel_mito_vars))$n_samples_conf_detected,
                     strand_correlation=rowMeans(cbind(data.frame(rowData(dx_mito_vars))$strand_correlation,data.frame(rowData(rel_mito_vars))$strand_correlation)),
                     vmr=rowMeans(cbind(data.frame(rowData(dx_mito_vars))$vmr,data.frame(rowData(rel_mito_vars))$vmr)),
                     mean_coverage=rowMeans(cbind(data.frame(rowData(dx_mito_vars))$mean_coverage,data.frame(rowData(rel_mito_vars))$mean_coverage)))
  rownames(mut_df)<-rownames(dx_mito_vars)
  
  # filter out low quality variants
  # filtered_variants<-mut_df %>% 
  #   filter(mut_df$n_cells_conf_detected >= 3 & 
  #            mut_df$strand_correlation > 0.65 & 
  #            log10(mut_df$vmr) > -3 & 
  #            mut_df$mean_coverage >= 10)
  # 
  # message(paste("Number of filtered variants:", nrow(filtered_variants)))
  
  dx_allele_freq<-assays(dx_mito_vars)[[1]] %>% as.matrix()
  rel_allele_freq<-assays(rel_mito_vars)[[1]] %>% as.matrix()
  
  # if(ncol(dx_allele_freq)<ncol(rel_allele_freq)){
  #   rel_allele_freq<-rel_allele_freq[,sample(colnames(rel_allele_freq), ncol(dx_allele_freq))]
  # }else{
  #   dx_allele_freq<-dx_allele_freq[,sample(colnames(dx_allele_freq), ncol(rel_allele_freq))]
  # }
  
  allele_freq<-cbind(dx_allele_freq,rel_allele_freq)
  # allele_freq<-allele_freq[rownames(filtered_variants),]
  
  allele_freq<-allele_freq[intersect(rownames(allele_freq),manual_mutations[[p]]),]
  
  allele_freq_binary<-allele_freq
  allele_freq_binary[allele_freq_binary>0]<-1
  
  # print(head(rev(sort(rowSums(allele_freq_binary))), 20))
  # hist(rowSums(allele_freq_binary), main = p)
  
  allele_freq_binary<-rbind(allele_freq_binary,rep(1,ncol(allele_freq_binary)))
  
  cl <- seuratSNN(sqrt(t(allele_freq_binary)), manual_resolution[[p]], 10)
  cl <- str_pad(cl, 2, pad = "0")
  table(cl)
  
  allele_freq_binary<-allele_freq_binary[rownames(allele_freq_binary)!="",]
  
  #cluster_ref<-data.frame(cluster=cl)
  #rownames(cluster_ref)<-colnames(allele_freq)
  #allele_freq<-reorder_dataframe(df=allele_freq, reference = cluster_ref, sample_order = sort(unique(cl)))
  
  
  # calculate avg heteroplasmy per cluster
  avg_heteroplasmy<-data.frame(matrix(nrow=nrow(allele_freq),ncol=0))
  for(c in sort(unique(cl))){
    message(c)
    x<-colnames(avg_heteroplasmy)
    avg_heteroplasmy<-cbind(avg_heteroplasmy,rowMeans(allele_freq_binary[,cl==c]))
    colnames(avg_heteroplasmy)<-c(x,c)
  }
  
  #var_to_keep<-rownames(avg_heteroplasmy)[rowMax(as.matrix(avg_heteroplasmy))>0.01]
  #avg_heteroplasmy<-avg_heteroplasmy[var_to_keep,]
  #allele_freq<-allele_freq[var_to_keep,]
  
  assignments<-colnames(avg_heteroplasmy)[max.col(avg_heteroplasmy)]
  feature_order<-c()
  for(c in colnames(avg_heteroplasmy)){
    regions<-rownames(avg_heteroplasmy)[assignments==c]
    feature_order<-c(feature_order,regions[rev(order(avg_heteroplasmy[regions,c]))])
  }
  allele_freq<-allele_freq[intersect(manual_mutations[[p]], rownames(allele_freq)), order(cl)]
  
  saveRDS(allele_freq, paste0("outputs/scATAC/mito_vars_data/mito_mutations/",p,"_mito_mutations_filtered.rds"))
  
  # Make annotations
  if(T){
    annotations<-data.frame(cluster=sort(cl),
                            timepoint=strsplit(colnames(allele_freq),"_") %>% sapply('[',1))
    rownames(annotations)<-colnames(allele_freq)
    annotations$cluster<-as.factor(annotations$cluster)
    annotations$timepoint<-as.factor(annotations$timepoint)
    #cluster_colors<-ggpubr::get_palette("ucscgb",length(unique(annotations$cluster)))
    cluster_colors<-ArchRPalettes$stallion[1:length(unique(annotations$cluster))]
    timepoint_colors<-c("darkblue","darkred")
    names(cluster_colors)<-levels(annotations$cluster)
    names(timepoint_colors)<-levels(annotations$timepoint)
    anno_colors<-list(cluster=cluster_colors,
                      timepoint=timepoint_colors)
  }
  
  heatmap_colors<-colorRampPalette(c("white", "#faee05", "#d40000"))
  
  h<-Heatmap(allele_freq,
             border = T,
             cluster_columns = F,
             cluster_rows = F,
             col=circlize::colorRamp2(seq(0.01,.05,length.out = 10), heatmap_colors(10)),
             show_row_names=T,
             show_column_names=F,
             heatmap_legend_param = list(title = "Heteroplasmy"),
             column_title=p,
             top_annotation=HeatmapAnnotation(df=annotations, col = anno_colors,  annotation_name_side = "left"),
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 8), 
             use_raster = T)
  print(h)
  
  pdf(paste0("outputs/scATAC_plots/", p, "_mito_heteroplasmy_heatmap.pdf"), width=6, height=3)
  print(h)
  dev.off()
  
  cluster_counts<-as.data.frame.matrix(table(annotations))
  cluster_counts$relapse_fraction<-cluster_counts$REL/rowSums(cluster_counts)
  cluster_counts$cluster<-rownames(cluster_counts)
  
  plt<-ggplot(cluster_counts, aes(x=cluster,y=relapse_fraction)) +
    geom_bar(stat = "identity") +
    theme_classic()
  
  #print(plt)
  
  saveRDS(tibble::rownames_to_column(annotations, "ID"), paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  
  allele_freq_melt<-as.data.frame(reshape2::melt(as.matrix(allele_freq)))
  colnames(allele_freq_melt)<-c("Var","ID","heteroplasmy")
  allele_freq_melt$cluster<-annotations$cluster[match(allele_freq_melt$ID,rownames(annotations))]
  
  plt_list<-list()
  for(v in unique(allele_freq_melt$Var)){
    plt<-ggplot(allele_freq_melt[allele_freq_melt$Var==v&allele_freq_melt$cluster%ni%as.character(11:100),], aes(x=cluster, y=heteroplasmy, fill=cluster)) +
      geom_violin() +
      theme_classic() +
      scale_fill_manual(values = as.vector(ArchRPalettes$stallion[1:11])) +
      ggtitle(paste(p,v)) +
      theme(legend.position = "none")
    plt_list[[v]]<-plt
  }
  pdf(paste0("outputs/scATAC_plots/", p,"_cluster_heteroplasmy_violin.pdf"), width=10, height=10)
  print(ggarrange(plotlist = plt_list))
  dev.off()
  
}


