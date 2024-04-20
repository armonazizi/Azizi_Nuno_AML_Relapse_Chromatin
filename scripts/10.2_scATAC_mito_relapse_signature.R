# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Get relapse signature for each mito clone
#
# Workflow for this script is as follows:
# 1. score each cell based on relapse gene accessibility signature
# 2. Plot violin of clones and relapse signature
# 3. 

library(ArchR)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")


# load and plot seaprate archr projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(ArchRProjects)<-unique(scATAC_ref$patient)

# do plotting
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
  
  mito_clone_data<-readRDS(paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  
  clones_to_keep<-table(mito_clone_data$cluster)
  clones_to_keep<-names(clones_to_keep)[clones_to_keep>50]
  
  mito_clone_data<-mito_clone_data[mito_clone_data$cluster%in%clones_to_keep,]
  
  ArchRProjects[[p]]@cellColData$timepoint<-rep("DX",nrow(ArchRProjects[[p]]@cellColData))
  ArchRProjects[[p]]@cellColData$timepoint[nchar(ArchRProjects[[p]]@cellColData$Sample)==6]<-"REL"
  ArchRProjects[[p]]@cellColData$ID<-paste(ArchRProjects[[p]]@cellColData$timepoint,strsplit(rownames(ArchRProjects[[p]]@cellColData),"#") %>% sapply(.,'[',2),sep="_")
  
  mito_clone_data$rel_sig<-ArchRProjects[[p]]@cellColData$rel_sig[match(mito_clone_data$ID,ArchRProjects[[p]]@cellColData$ID)]
  
  mito_clone_data$clone_timepoint<-paste(mito_clone_data$cluster,mito_clone_data$timepoint,sep="_")
  
  colors_transparent<-col2rgb(ArchRPalettes$stallion)
  colors_transparent<-rgb(colors_transparent["red",],colors_transparent["green",],colors_transparent["blue",],maxColorValue = 255, alpha = .4*255)
  new_pal<-c(rbind(colors_transparent,ArchRPalettes$stallion))
  names(new_pal)<-c(rbind(paste0(str_pad(0:19, 2, pad = "0"),"_DX"),paste0(str_pad(0:19, 2, pad = "0"),"_REL")))
  
  #plot violin
  p1 <- ggplot(mito_clone_data, aes(x=clone_timepoint, y=rel_sig, color=clone_timepoint, fill=clone_timepoint)) +
    #geom_boxplot(color="black",outlier.alpha = 0, show.legend = F) +
    geom_violin(color="black") +
    theme_classic() +
    scale_color_manual(values=new_pal) +
    scale_fill_manual(values=new_pal) +
    ggrastr::rasterize(geom_jitter(width = 0.1, size=0.2, show.legend = F, color="black"), dpi=300) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black",size=10)) +
    xlab("") +
    ylab("Relapse Signature Score")
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_mito_clone_violin_relapse_signature_peaks_enrichment.pdf"), width = 4, height = 2.5)
  print(p1)
  dev.off()
  
  #plot boxplot
  p1 <- ggplot(mito_clone_data, aes(x=clone_timepoint, y=rel_sig, color=clone_timepoint, fill=clone_timepoint)) +
    geom_boxplot(color="black",outlier.alpha = 0) +
    #geom_violin(color="black") +
    theme_classic() +
    scale_color_manual(values=new_pal) +
    scale_fill_manual(values=new_pal) +
    ggrastr::rasterize(geom_jitter(width = 0.1, size=0.2, show.legend = F, color="black"), dpi=300) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black",size=10)) +
    xlab("") +
    ylab("Relapse Signature Score")
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_mito_clone_boxplot_relapse_signature_peaks_enrichment.pdf"), width = 4.5, height = 2.5)
  print(p1)
  dev.off()
  
  
  #timepoint boxplot
  p1 <- ggplot(mito_clone_data, aes(x=timepoint, y=rel_sig, color=timepoint, fill=timepoint)) +
    geom_boxplot(color="black",outlier.alpha = 0) +
    #geom_violin(color="black") +
    theme_classic() +
    scale_color_manual(values=c("darkblue","darkred")) +
    scale_fill_manual(values=c("darkblue","darkred")) +
    ggrastr::rasterize(geom_jitter(width = 0.1, size=0.1, show.legend = F, color="black"), dpi=300) +
    theme(axis.text.x = element_text(color="black",size=10)) +
    xlab("") +
    ylab("Relapse Signature Score")
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_timepoint_boxplot_relapse_signature_peaks_enrichment.pdf"), width = 3, height = 2.5)
  print(p1)
  dev.off()
  
  
  #timepoint density plot
  p1 <- ggplot(mito_clone_data, aes(x=rel_sig, y=timepoint, color=timepoint, fill=timepoint)) +
    geom_density_ridges(aes(fill = timepoint)) +
    theme_classic() +
    scale_color_manual(values=c("darkblue","darkred")) +
    scale_fill_manual(values=c("darkblue","darkred")) +
    theme(axis.text.x = element_text(color="black",size=10)) +
    ylab("") +
    xlab("Relapse Signature Score")
  
  #visualize two plots above side-by-side
  pdf(paste0("outputs/scATAC_plots/",p,"_scATAC_timepoint_distribution_relapse_signature_peaks_enrichment.pdf"), width = 4, height = 1.5)
  print(p1)
  dev.off()
  
}

