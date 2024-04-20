# Armon Azizi (aazizi@stanford.edu)
# Functions for Nuno/Azizi et al. for analyzing AML Relapse ATACseq data

# packages
if(!require(edgeR)){BiocManager::install("edgeR")}
if(!require(preprocessCore)){BiocManager::install("preprocessCore")}
if(!require(DESeq2)){BiocManager::install("DESeq2")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(tibble)){install.packages("tibble")}
if(!require(ComplexHeatmap)){BiocManager::install("ComplexHeatmap")}
if(!require(circlize)){install.packages("circlize")}
if(!require(viridis)){install.packages("viridis")}
if(!require(fgsea)){BiocManager::install("fgsea")}
if(!require(data.table)){install.packages("data.table")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(plyr)){install.packages("plyr")}
if(!require(ggridges)){install.packages("ggridges")}

library(edgeR)
library(preprocessCore)
library(DESeq2)
library(ggplot2)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(fgsea)
library(data.table)
library(openxlsx)
library(plyr)
library(ggridges)
library(stringr)

# Normalize count matrix using edgeR cpm() and preprocessCore normalize.quantiles
# @param count_matrix dataframe or matrix with samples as columns and genes/peaks as rows
# @param log should the data be transformed to log space when performing cpm?
# @param quantile should the data be quantile normalized
# @return normalized dataframe
# Same normalization as used in corces et al. Science
normalize_count_matrix<-function(count_matrix, log=TRUE, quantile=TRUE){
  count_matrix<-cpm(count_matrix, log = log, prior.count=5)
  if(quantile){
    temp_cm<-as.data.frame(normalize.quantiles(count_matrix))
    colnames(temp_cm)<-colnames(count_matrix)
    rownames(temp_cm)<-rownames(count_matrix) 
  }else{
    temp_cm<-count_matrix
  }
  return(temp_cm)
}

# function to collapse replicates
collapse_count_matrix<-function(count_matrix, sample_ref, id="sample"){
  result<-count_matrix[,c()]
  for(s in unique(sample_ref[,id])){
    result[,s]<-rowMeans(cbind(count_matrix[,sample_ref$name[sample_ref[,id]==s&sample_ref$replicate==1]],
                               count_matrix[,sample_ref$name[sample_ref[,id]==s&sample_ref$replicate==2]]))
  }
  return(result)
}

# Clustering function based in seurat
seuratSNN <- function(matSVD, resolution = 1, k.param = 10){ 
  set.seed(1)
  rownames(matSVD) <- make.unique(rownames(matSVD))
  obj <- FindNeighbors(matSVD, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

# Function to reorder a dataframe for clean heatmap plotting
reorder_dataframe<-function(df,reference,sample_order){
  sample_order<-as.factor(sample_order)
  colnames(reference)<-c("CellType")
  reference$CellType<-as.factor(reference$CellType)
  #reference$CellType<-reorder(reference$CellType, new.order=sample_order)
  so<-rownames(reference[order(reference$CellType),,drop=FALSE])
  
  collapsed_counts<-data.frame(matrix(nrow=nrow(df),ncol=0))
  
  for(celltype in sample_order){
    message(celltype)
    x<-colnames(collapsed_counts)
    collapsed_counts<-cbind(collapsed_counts,rowMeans(df[,rownames(reference)[reference$CellType==celltype]]))
    colnames(collapsed_counts)<-c(x,celltype)
  }
  assignments<-colnames(collapsed_counts)[max.col(collapsed_counts)]
  feature_order<-c()
  for(celltype in colnames(collapsed_counts)){
    regions<-rownames(collapsed_counts)[assignments==celltype]
    feature_order<-c(feature_order,regions[rev(order(collapsed_counts[regions,celltype]))])
  }
  return(df[feature_order,so])
}

# --------------------------------------------------
# get variable genes
# --------------------------------------------------
# m: normalized matrix
get_variable_gene<-function(m) {
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,unique(quantile(mean,seq(0.1,1,0.05)),Inf))))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}


reorder_dataframe<-function(df,reference,sample_order){
  sample_order<-as.factor(sample_order)
  colnames(reference)<-c("CellType")
  reference$CellType<-as.factor(reference$CellType)
  reference$CellType<-factor(reference$CellType, levels =sample_order)
  so<-rownames(reference[order(reference$CellType),,drop=FALSE])
  
  collapsed_counts<-data.frame(matrix(nrow=nrow(df),ncol=0))
  
  for(celltype in sample_order){
    message(celltype)
    x<-colnames(collapsed_counts)
    collapsed_counts<-cbind(collapsed_counts,rowMeans(df[,rownames(reference)[reference$CellType==celltype],drop=F]))
    colnames(collapsed_counts)<-c(x,celltype)
  }
  assignments<-colnames(collapsed_counts)[max.col(collapsed_counts)]
  feature_order<-c()
  for(celltype in colnames(collapsed_counts)){
    regions<-rownames(collapsed_counts)[assignments==celltype]
    feature_order<-c(feature_order,regions[rev(order(collapsed_counts[regions,celltype]))])
  }
  return(df[feature_order,so])
}

get_feature_assignments<-function(df,reference,sample_order){
  sample_order<-as.factor(sample_order)
  colnames(reference)<-c("CellType")
  reference$CellType<-as.factor(reference$CellType)
  reference$CellType<-factor(reference$CellType,levels=sample_order)
  so<-rownames(reference[order(reference$CellType),,drop=FALSE])
  
  
  collapsed_counts<-data.frame(matrix(nrow=nrow(df),ncol=0))
  
  for(celltype in sample_order){
    message(celltype)
    x<-colnames(collapsed_counts)
    collapsed_counts<-cbind(collapsed_counts,rowMeans(df[,rownames(reference)[reference$CellType==celltype],drop=F]))
    colnames(collapsed_counts)<-c(x,celltype)
  }
  assignments<-colnames(collapsed_counts)[max.col(collapsed_counts)]
  
  return(data.frame(feature=rownames(df),assignment=assignments))
}

