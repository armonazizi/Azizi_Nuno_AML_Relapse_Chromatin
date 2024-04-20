# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Plotting and analysis of scATAC data from DX and REL AML
#
# Workflow for this script is as follows:
# 1. Load pregenerated ArchR projects
# 2. Plot umaps for all samples both combined and independently
# 3. 
# 4. 

library(ArchR)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# Load combined ArchR project
ArchR_comb<-loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/all_samples_combined"))

#plot UMAP
p1 <- plotEmbedding(ArchRProj = ArchR_comb, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_postfilter", 
                    labelMeans=F,
                    rastr = T, dpi=250)

#plot UMAP colored by Cluster instead of Sample
p2 <- plotEmbedding(ArchRProj = ArchR_comb, 
                    colorBy = "cellColData", 
                    name = "timepoint", 
                    pal = c("#272E6A","#D51F26"), 
                    embedding = "UMAP_postfilter", 
                    labelMeans=F,
                    rastr = T, dpi=250)

#visualize two plots above side-by-side
# print(ggAlignPlots(p7, p8, type = "h", draw = T))

pdf("outputs/scATAC_plots/scATAC_combined_samples_umap.pdf", width = 5, height = 5)
p1
dev.off()

pdf("outputs/scATAC_plots/scATAC_combined_samples_umap_timepoint_labels.pdf", width = 5, height = 5)
p2
dev.off()


colors_transparent<-col2rgb(ArchRPalettes$stallion)
colors_transparent<-rgb(colors_transparent["red",],colors_transparent["green",],colors_transparent["blue",],maxColorValue = 255, alpha = .1*255)
new_pal<-c(rbind(colors_transparent,ArchRPalettes$stallion))[1:length(unique(ArchR_comb@cellColData$Sample))]
new_pal<-get_palette("Paired",8)
new_pal<-c("#fc9a9a","#cc0000","#a8ecff","#0a00c7","#e88cfa","#6202bd","#bdff9c","#0a7d00")
names(new_pal)<-unique(ArchR_comb@cellColData$Sample)

#plot UMAP new colors
p1 <- plotEmbedding(ArchRProj = ArchR_comb, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_postfilter", 
                    randomize = T,
                    pal = new_pal, 
                    labelMeans=F,
                    rastr = T, dpi=250)
print(p1)

pdf("outputs/scATAC_plots/scATAC_combined_samples_umap_paired_colors.pdf", width = 5, height = 5)
p1
dev.off()


legend <- ggplot(data.frame(Sample=names(new_pal),b=factor(1:8)), aes(Sample, fill = Sample)) + 
  geom_bar() +
  scale_fill_manual(values=new_pal)

pdf("outputs/scATAC_plots/scATAC_combined_samples_umap_paired_colors_legend.pdf", width = 5, height = 5)
legend
dev.off()



# load and plot seaprate archr projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(ArchRProjects)<-unique(scATAC_ref$patient)

# do umaps
for(p in unique(scATAC_ref$patient)){
  message(p)
  
  #plot UMAP
  p1 <- plotEmbedding(ArchRProj = ArchRProjects[[p]], 
                      colorBy = "cellColData", 
                      name = "Sample", 
                      pal = c("#272E6A","#D51F26"), 
                      embedding = "UMAP_postfilter",
                      rastr = T, dpi=250)
  #plot UMAP colored by Cluster instead of Sample
  p2 <- plotEmbedding(ArchRProj = ArchRProjects[[p]], 
                      colorBy = "cellColData", 
                      name = "Clusters", 
                      embedding = "UMAP_postfilter",
                      rastr = T, dpi=250)
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_umap_timepoint_labels.pdf"), width = 5, height = 5)
  print(p1)
  dev.off()
  
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_umap_cluster_labels.pdf"), width = 5, height = 5)
  print(p2)
  dev.off()
}
