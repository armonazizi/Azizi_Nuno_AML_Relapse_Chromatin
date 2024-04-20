# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Do analyses of gene and TF accessibility across timepoints and pseudotime
#
# Workflow for this script is as follows:
# 1. 
# 2. 
# 3. 
# 4. 


library(ArchR)
library(openxlsx)
library(Seurat)
library(ggpubr)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
scATAC_ref$timepoint<-rep("DX", nrow(scATAC_ref))
scATAC_ref$timepoint[nchar(scATAC_ref$ID)==6]<-"REL"
scATAC_ref$name<-paste(scATAC_ref$patient, scATAC_ref$timepoint, sep="_")

ArchRProjects<-lapply(unique(scATAC_ref$name), function(x) loadArchRProject(paste0("~/Bioinformatics/AML_relapse_project/analysis_clean/external_input_files/ArchR_Projects/AML_ArchR_Projects_hg38_mito/",x)))
names(ArchRProjects)<-unique(scATAC_ref$name)

p<-"SU360"

# label LSCs
if(file.exists(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))){
  lsc_classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_LSC_classifications/",p,"_scATAC_lsc_classifications.rds"))
  
  ArchRProjects[[paste0(p,"_DX")]]@cellColData$CellType<-lsc_classifications$classification[match(rownames(ArchRProjects[[paste0(p,"_DX")]]@cellColData),rownames(lsc_classifications))]
  ArchRProjects[[paste0(p,"_REL")]]@cellColData$CellType<-lsc_classifications$classification[match(rownames(ArchRProjects[[paste0(p,"_REL")]]@cellColData),rownames(lsc_classifications))]
  
  plt1<-plotEmbedding(ArchRProj = ArchRProjects[[paste0(p,"_DX")]], 
                      colorBy = "cellColData", 
                      name = "CellType",
                      embedding = "UMAP", 
                      labelMeans=F)
  
  plt2<-plotEmbedding(ArchRProj = ArchRProjects[[paste0(p,"_REL")]], 
                      colorBy = "cellColData", 
                      name = "CellType",
                      embedding = "UMAP", 
                      labelMeans=F)
  pdf(paste0("outputs/scATAC_plots/pseudotime/",p,"_scATAC_umap_LSC_classification_labels.pdf"), width = 8, height = 5)
  print(ggarrange(plt1,plt2, nrow=1))
  dev.off()
}


# Label projected celltypes
classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
ArchRProjects[[paste0(p,"_DX")]]@cellColData$projection_classification<-classifications$AML_classification_original_labels[match(rownames(ArchRProjects[[paste0(p,"_DX")]]@cellColData),rownames(classifications))]
ArchRProjects[[paste0(p,"_REL")]]@cellColData$projection_classification<-classifications$AML_classification_original_labels[match(rownames(ArchRProjects[[paste0(p,"_REL")]]@cellColData),rownames(classifications))]

plt1<-plotEmbedding(ArchRProj = ArchRProjects[[paste0(p,"_DX")]], 
                   colorBy = "cellColData", 
                   name = "projection_classification",
                   embedding = "UMAP", size=0.5,
                   labelMeans=T)

plt2<-plotEmbedding(ArchRProj = ArchRProjects[[paste0(p,"_REL")]], 
                    colorBy = "cellColData", 
                    name = "projection_classification",
                    embedding = "UMAP", size=0.5,
                    labelMeans=T)
pdf(paste0("outputs/scATAC_plots/pseudotime/",p,"_scATAC_umap_projection_classification_labels.pdf"), width = 8, height = 5)
print(ggarrange(plt1,plt2, nrow=1))
dev.off()

# add trajectory
dx_trajectory <- c("05_CMP.LMPP", "GMP-08_GMP.Neut", "07_GMP")
rel_trajectory <-  c("05_CMP.LMPP", "GMP-08_GMP.Neut", "07_GMP")

ArchRProjects[[paste0(p,"_DX")]] <- addTrajectory(
  ArchRProj = ArchRProjects[[paste0(p,"_DX")]], 
  name = "CMP_LMPP_to_GMP", 
  groupBy = "projection_classification",
  trajectory = dx_trajectory, 
  embedding = "UMAP",spar = 1.5,
  force = TRUE
)

ArchRProjects[[paste0(p,"_REL")]] <- addTrajectory(
  ArchRProj = ArchRProjects[[paste0(p,"_REL")]], 
  name = "CMP_LMPP_to_GMP", 
  groupBy = "projection_classification",
  trajectory = rel_trajectory,
  embedding = "UMAP", spar = 1.5,
  force = TRUE
)

plt1<-plotTrajectory(ArchRProjects[[paste0(p,"_DX")]], trajectory = "CMP_LMPP_to_GMP", colorBy = "cellColData", name = "CMP_LMPP_to_GMP")
plt2<-plotTrajectory(ArchRProjects[[paste0(p,"_REL")]], trajectory = "CMP_LMPP_to_GMP", colorBy = "cellColData", name = "CMP_LMPP_to_GMP")

pdf(paste0("outputs/scATAC_plots/pseudotime/",p,"_scATAC_pseudotime_trajectory.pdf"), width = 8, height = 5)
print(ggarrange(plt1[[1]],plt2[[1]]))
dev.off()



# ArchRProjects[[paste0(p,"_DX")]] <- addImputeWeights(ArchRProjects[[paste0(p,"_DX")]])
# ArchRProjects[[paste0(p,"_REL")]] <- addImputeWeights(ArchRProjects[[paste0(p,"_REL")]])
# 
# 
# plt1 <- plotTrajectory(ArchRProjects[[paste0(p,"_DX")]], trajectory = "CMP_LMPP_to_GMP", colorBy = "GeneScoreMatrix", name = "GATA2", continuousSet = "horizonExtra")
# plt2 <- plotTrajectory(ArchRProjects[[paste0(p,"_REL")]], trajectory = "CMP_LMPP_to_GMP", colorBy = "GeneScoreMatrix", name = "GATA2", continuousSet = "horizonExtra")
# print(ggarrange(plt1[[1]],plt2[[1]],plt1[[2]],plt2[[2]]), nrow=2)
# 
# 
# plt_data  <- getTrajectory(ArchRProj = archr_sub, name = "CMP_LMPP_to_GMP", useMatrix = "GeneScoreMatrix", log2Norm = T)
# labels<-grep(":CD[0-9]",rownames(plt_data), value = T)
# plt <- plotTrajectoryHeatmap(plt_data,  pal = paletteContinuous(set = "horizonExtra"), labelMarkers = labels, labelTop=0)
# print(plt)


trajMM  <- getTrajectory(ArchRProj = ArchRProjects[[paste0(p,"_DX")]], name = "CMP_LMPP_to_GMP", useMatrix = "MotifMatrix", log2Norm = FALSE, scaleTo = NULL)
mtx<-plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), returnMatrix = T)
plt1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
#print(plt)

trajMM  <- getTrajectory(ArchRProj = ArchRProjects[[paste0(p,"_REL")]], name = "CMP_LMPP_to_GMP", useMatrix = "MotifMatrix", log2Norm = FALSE, scaleTo = NULL)
plt2 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), useSeqnames = rownames(mtx), rowOrder = rownames(mtx))
#print(plt)

pdf(paste0("outputs/scATAC_plots/pseudotime/",p,"_scATAC_pseudotime_TF_heatmap.pdf"), width = 10, height = 7)
print(plt1+plt2)
dev.off()



#plt<-plotTrajectory(ArchRProjects[[paste0(p,"_DX")]], trajectory = "CMP_LMPP_to_GMP", colorBy = "MotifMatrix", name = "z:CEBPA_155", continuousSet = "blueYellow")
