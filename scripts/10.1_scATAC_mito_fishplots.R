# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# Generate fishplots of the mitochondrial clonal evolution hierarchy
#
#
# Workflow for this script is as follows:
# 1. 
# 2. 
# 3. 

library(fishplot)
library(stringr)
library(dplyr)
library(ArchR)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
scATAC_ref<-readRDS("inputs/scATAC/sample_reference.rds")
scATAC_ref$timepoint<-rep("DX", nrow(scATAC_ref))
scATAC_ref$timepoint[nchar(scATAC_ref$ID)==6]<-"REL"

# make fishplot for each patient

parents_list<-list(SU142 = c(0,1,1,2),
                   SU360 = c(5,1,1,3,0),
                   SU484 = c(0,1,1,2,1,2,1),
                   SU892 = c(3,0,2,2))
origin_clones<-list(SU142 = c(1),
                    SU360 = c(5),
                    SU484 = c(1),
                    SU892 = c(2))
cutoff<-list(SU142 = 1,
             SU360 = 2,
             SU332 = 1,
             SU484 = 6,
             SU892 = 2)

for(p in unique(scATAC_ref$patient)){
  message(p)
  mito_clusters<-readRDS(paste0("outputs/scATAC/mito_vars_data/mito_clones/",p,"_mito_clones.rds"))
  mito_clusters$ID<-gsub("DX_",paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="DX"],"#"), mito_clusters$ID)
  mito_clusters$ID<-gsub("REL_",paste0(scATAC_ref$ID[scATAC_ref$patient==p&scATAC_ref$timepoint=="REL"],"#"), mito_clusters$ID)
  
  cluster_fractions<-as.data.frame.matrix(table(mito_clusters[c("cluster","timepoint")]))
  cluster_fractions<-sweep(cluster_fractions, 2, colSums(cluster_fractions),"/")
  cluster_fractions<-cluster_fractions*100
  cluster_fractions<-cluster_fractions[rowSums(cluster_fractions)>cutoff[[p]],]
  
  parents = parents_list[[p]]
  
  # for(i in 1:length(parents)){
  #   if(parents[i]!=0){
  #     cluster_fractions[parents[i],]<-cluster_fractions[parents[i],]+cluster_fractions[i,]
  #   }
  # }
  
  add_parents<-function(fractions, parents, idx){
    children<-which(parents==idx)
    message(children)
    if(length(children)==0){return(fractions)}
    
    for(c in children){
      fractions<-add_parents(fractions, parents, c)
      fractions[idx,]<-fractions[idx,]+fractions[c,]
    }
    
    return(fractions)
  }
  
  origins<-origin_clones[[p]]
  for(o in origins){
    cluster_fractions<-add_parents(cluster_fractions, parents, o)
  }
  
  #provide a list of timepoints to plot
  #You may need to add interpolated points to end up with the desired
  #visualization. Example here was actually sampled at days 0 and 150
  timepoints=c(0,50,100,150)      
  
  #provide a matrix with the fraction of each population
  #present at each timepoint
  
  frac.table<-cluster_fractions
  frac.table<-cbind(frac.table[,1],cbind(frac.table[,1]/100,cbind(frac.table[,2]/100,frac.table[,2])))
  frac.table<-frac.table*0.99
  
  #provide a vector listing each clone's parent
  #(0 indicates no parent)
  parents = parents_list[[p]]
  
  #create a fish object
  fish = createFishObject(frac.table,parents,timepoints=timepoints, col = ArchRPalettes$stallion[1:nrow(frac.table)])
  
  #calculate the layout of the drawing
  fish = layoutClones(fish)
  
  #draw the plot, using the splining method (recommended)
  #and providing both timepoints to label and a plot title
  pdf(paste0("outputs/scATAC_plots/",p,"_mito_fishplot.pdf"), width=5, height=3.5)
  print(fishPlot(fish,shape="spline",title.btm=p,
           cex.title=0.5, vlines=c(0,150), 
           vlab=c("DX","REL"), 
           bg.type = "solid", bg.col = "grey90"))
  dev.off()
}
