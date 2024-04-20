# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 5
# Projection of mtscATAC to healthy manifold
#

library(ArchR)
library(openxlsx)
library(Seurat)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference_hg19.rds")
scATAC_ref$timepoint<-rep("DX", nrow(scATAC_ref))
scATAC_ref$timepoint[nchar(scATAC_ref$ID)==6]<-"REL"

AML_ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg19_mito/",x)))
names(AML_ArchRProjects)<-unique(scATAC_ref$patient)


for(p in unique(scATAC_ref$patient)){
  LSI_projection<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  
  mito_clones<-readRDS(paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  
  mito_clones$ID[mito_clones$timepoint=="DX"]<-paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="DX"], "#",
                                                    sapply(strsplit(mito_clones$ID[mito_clones$timepoint=="DX"],"_"),'[',2))
  
  mito_clones$ID[mito_clones$timepoint=="REL"]<-paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="REL"], "#",
                                                      sapply(strsplit(mito_clones$ID[mito_clones$timepoint=="REL"],"_"),'[',2))
  
  LSI_projection$clone<-as.character(mito_clones$cluster)[match(rownames(LSI_projection),mito_clones$ID)]
  LSI_projection$clone[LSI_projection$Sample=="reference"]<-"reference"
  LSI_projection<-LSI_projection[complete.cases(LSI_projection),]
  
  p1<-ggPoint(LSI_projection$X1, LSI_projection$X2, color = LSI_projection$clone, 
              pal = c(ArchRPalettes$stallion[1:(length(unique(LSI_projection$clone))-1)],"lightgrey"), 
              #colorOrder = c(names(celltype_pal[names(celltype_pal)%in%unique(projectionDF$AML_classification)]),"reference"),
              labelMeans = F, 
              xlabel = "UMAP Dimension 1",
              ylabel = "UMAP Dimension 2", 
              size=0.1,
              baseSize = 6, 
              labelSize = 1.5,
              rastr = T, dpi=250)
  pdf(paste0("outputs/scATAC_plots/scATAC_mito_projection/",p,"_mito_clone_projection.pdf"), width = 4, height = 4, useDingbats = FALSE)
  print(p1)
  dev.off()
}




files<-list.files("inputs/scATAC/mito_vars_data/") %>% grep("rds", ., value = T)

file_reference<-data.frame(patient=substr(files,0,5),
                           sample=sapply(strsplit(files, "-|_"),'[',1),
                           file=files)
file_reference$timepoint<-rep("DX",nrow(file_reference))
file_reference$timepoint[nchar(file_reference$sample)==6]<-"REL"


manual_mutations<-list(
  SU142=c("16215A>G", "933G>A", "13369T>C", "4429G>A",  "3496G>A", "1107A>G"),
  SU332=c("13129C>T", "11719G>A", "8251G>A", "11928A>G", "3921C>A", "1719G>A", "1703C>T", "14581T>C", "4960C>T", "16131T>C", "9310T>C"),
  SU360=c("2471G>A", "6810G>A", "11838T>A", "2805A>T", "7775G>A", "6456G>A", "12871G>A", "1040T>C", "7849C>T", "2908T>C"),
  SU484=c("12476G>A", "3244G>A", "7069T>C", "15713T>C", "3199T>C", "6472T>C"),
  SU892=c("15639T>G", "16189T>C"),
  SU926=c("9210A>G", "12414T>C", "7061A>G", "11955A>G", "523A>C")
)


for(p in unique(scATAC_ref$patient)){
  LSI_projection<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  
  message(p)
  
  dx_mito_vars<-readRDS(paste0("inputs/scATAC/mito_vars_data/",file_reference$file[file_reference$patient==p&file_reference$timepoint=="DX"]))
  rel_mito_vars<-readRDS(paste0("inputs/scATAC/mito_vars_data/",file_reference$file[file_reference$patient==p&file_reference$timepoint=="REL"]))
  
  colnames(dx_mito_vars)<-paste0("DX_",colnames(dx_mito_vars))
  colnames(rel_mito_vars)<-paste0("REL_",colnames(rel_mito_vars))
  
  dx_allele_freq<-assays(dx_mito_vars)[[1]] %>% as.matrix()
  rel_allele_freq<-assays(rel_mito_vars)[[1]] %>% as.matrix()
  
  allele_freq<-cbind(dx_allele_freq,rel_allele_freq)
  
  allele_freq<-allele_freq[intersect(rownames(allele_freq),manual_mutations[[p]]),] %>% t()
  
  mito_clones<-readRDS(paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  
  rownames(mito_clones)<-mito_clones$ID
  
  mito_clones$ID[mito_clones$timepoint=="DX"]<-paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="DX"], "#",
                                                      sapply(strsplit(mito_clones$ID[mito_clones$timepoint=="DX"],"_"),'[',2))
  
  mito_clones$ID[mito_clones$timepoint=="REL"]<-paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="REL"], "#",
                                                       sapply(strsplit(mito_clones$ID[mito_clones$timepoint=="REL"],"_"),'[',2))
  
  rownames(allele_freq)<-mito_clones$ID[match(rownames(allele_freq),rownames(mito_clones))]
  allele_freq[allele_freq>0]<-1
  
  plt_list<-list()
  for(mut in colnames(allele_freq)){
    LSI_projection$mutation<-as.character(allele_freq[,mut])[match(rownames(LSI_projection),rownames(allele_freq))]
    LSI_projection$mutation[LSI_projection$Sample=="reference"]<-"reference"
    LSI_projection$mutation[LSI_projection$mutation=="0"]<-"WT"
    LSI_projection$mutation[is.na(LSI_projection$mutation)]<-"WT"
    LSI_projection$mutation[LSI_projection$mutation=="1"]<-"Mut"
    LSI_projection$mutation<-factor(LSI_projection$mutation, levels = c("reference","WT","Mut"))
    LSI_projection<-LSI_projection[order(LSI_projection$mutation),]
    LSI_projection$mutation<-as.character(LSI_projection$mutation)
    
    p1<-ggPoint(LSI_projection$X1, LSI_projection$X2, color = LSI_projection$mutation, 
                pal = c("red","grey80","grey30"), 
                #colorOrder = c(names(celltype_pal[names(celltype_pal)%in%unique(projectionDF$AML_classification)]),"reference"),
                labelMeans = F, 
                xlabel = "UMAP Dimension 1",
                ylabel = "UMAP Dimension 2", 
                size=0.05,
                baseSize = 6, 
                labelSize = 1.5, title = mut,
                rastr = T, dpi=250, randomize = F)
    plt_list[[mut]]<-p1
  }
  
  #ggarrange(plotlist = plt_list)
  
  pdf(paste0("outputs/scATAC_plots/scATAC_mito_projection/",p,"_mito_separate_mutation_projection.pdf"), width = 10, height = 10, useDingbats = FALSE)
  print(ggarrange(plotlist = plt_list))
  dev.off()
}



mutations_to_plot<-list(
  SU142=c("933G>A", "13369T>C", "4429G>A",  "3496G>A"),
  SU360=c("1040T>C"),
  SU484=c("7069T>C", "15713T>C"),
  SU892=c("16189T>C")
)


for(p in unique(scATAC_ref$patient)){
  LSI_projection<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  
  message(p)
  
  dx_mito_vars<-readRDS(paste0("inputs/scATAC/mito_vars_data/",file_reference$file[file_reference$patient==p&file_reference$timepoint=="DX"]))
  rel_mito_vars<-readRDS(paste0("inputs/scATAC/mito_vars_data/",file_reference$file[file_reference$patient==p&file_reference$timepoint=="REL"]))
  
  colnames(dx_mito_vars)<-paste0("DX_",colnames(dx_mito_vars))
  colnames(rel_mito_vars)<-paste0("REL_",colnames(rel_mito_vars))
  
  dx_allele_freq<-assays(dx_mito_vars)[[1]] %>% as.matrix()
  rel_allele_freq<-assays(rel_mito_vars)[[1]] %>% as.matrix()
  
  allele_freq<-cbind(dx_allele_freq,rel_allele_freq)
  
  allele_freq<-allele_freq[intersect(rownames(allele_freq),mutations_to_plot[[p]]),]
  allele_freq<-t(allele_freq)
  
  mito_clones<-readRDS(paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  
  rownames(mito_clones)<-mito_clones$ID
  
  mito_clones$ID[mito_clones$timepoint=="DX"]<-paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="DX"], "#",
                                                      sapply(strsplit(mito_clones$ID[mito_clones$timepoint=="DX"],"_"),'[',2))
  
  mito_clones$ID[mito_clones$timepoint=="REL"]<-paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="REL"], "#",
                                                       sapply(strsplit(mito_clones$ID[mito_clones$timepoint=="REL"],"_"),'[',2))
  
  rownames(allele_freq)<-mito_clones$ID[match(rownames(allele_freq),rownames(mito_clones))]
  allele_freq[allele_freq>0]<-1
  
  LSI_projection$mutation<-rep("reference", nrow(LSI_projection))
  
  for(mut in colnames(allele_freq)){
    LSI_projection$mutation[rownames(LSI_projection)%in%rownames(allele_freq)[allele_freq[,mut]==1]]<-mut
  }
  
  LSI_projection$mutation<-factor(LSI_projection$mutation, levels = c("reference",colnames(allele_freq)))
  LSI_projection<-LSI_projection[order(LSI_projection$mutation),]
  LSI_projection$mutation<-as.character(LSI_projection$mutation)
  
  p1<-ggPoint(LSI_projection$X1, LSI_projection$X2, color = LSI_projection$mutation, 
              pal = c(ArchRPalettes$kelly[1:ncol(allele_freq)],"grey80"), 
              #colorOrder = c(names(celltype_pal[names(celltype_pal)%in%unique(projectionDF$AML_classification)]),"reference"),
              labelMeans = F, 
              xlabel = "UMAP Dimension 1",
              ylabel = "UMAP Dimension 2", 
              size=0.05,
              baseSize = 6, 
              labelSize = 1.5, title=p,
              rastr = T, dpi=250, randomize = F)
  p1
  
  
  pdf(paste0("outputs/scATAC_plots/scATAC_mito_projection/",p,"_mito_mutation_projection.pdf"), width = 4, height = 4, useDingbats = FALSE)
  print(p1)
  dev.off()
  
  
  LSI_projection$timepoint<-rep("DX",nrow(LSI_projection))
  LSI_projection$timepoint[LSI_projection$Sample=="reference"]<-"reference"
  LSI_projection$timepoint[nchar(LSI_projection$Sample)==6]<-"REL"
  
  plt_input<-as.data.frame(table(LSI_projection[,c("timepoint","mutation")]))
  plt_input<-plt_input[plt_input$timepoint!="reference",]
  plt_input<-plt_input[plt_input$mutation!="reference",]
  
  pal<-c("darkred","darkblue")
  names(pal)<-c("DX","REL")
  
  plt<-ggplot(plt_input, aes(fill=timepoint, y=Freq, x=mutation)) + 
    geom_bar(position="fill", stat="identity", color="black", size=0.5) +
    scale_fill_manual(values = pal) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Mitochondrial Mutation") +
    ylab("Sample Frequency") +
    ggtitle(p)
  
  pdf(paste0("outputs/scATAC_plots/scATAC_mito_projection/",p,"_mito_mutation_sample_frequencies.pdf"), width = 3.5, height = 4, useDingbats = FALSE)
  print(plt)
  dev.off()
  
  
  LSI_projection$timepoint_mutation<-paste(LSI_projection$timepoint,LSI_projection$mutation,sep="_")
  LSI_projection$timepoint_mutation[LSI_projection$timepoint=="reference"]<-"reference"
  LSI_projection$timepoint_mutation[LSI_projection$mutation=="reference"]<-"reference"
  
  p1<-ggPoint(LSI_projection$X1, LSI_projection$X2, color = LSI_projection$timepoint_mutation, 
              pal = c(ArchRPalettes$kelly[1:(ncol(allele_freq)*2)],"grey80"), 
              #colorOrder = c(names(celltype_pal[names(celltype_pal)%in%unique(projectionDF$AML_classification)]),"reference"),
              labelMeans = F, 
              xlabel = "UMAP Dimension 1",
              ylabel = "UMAP Dimension 2", 
              size=0.05,
              baseSize = 6, 
              labelSize = 1.5, title=p,
              rastr = T, dpi=250, randomize = F)
  p1
  

}

