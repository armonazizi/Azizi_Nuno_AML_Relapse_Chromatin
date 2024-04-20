# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 3
# Projection of scATAC to healthy manifold
#

library(ArchR)
library(openxlsx)
library(Seurat)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# load ArchR projects
heme_archr<-loadArchRProject("external_input_files/ArchR_Projects/healthy_hematopoiesis_ArchR")

scATAC_ref<-readRDS("inputs/scATAC/sample_reference_hg19.rds")

AML_ArchRProjects<-lapply(unique(scATAC_ref$patient), function(x) loadArchRProject(paste0("external_input_files/ArchR_Projects/AML_ArchR_Projects_hg19_mito/",x)))
names(AML_ArchRProjects)<-unique(scATAC_ref$patient)


# Add peak matrices for each AML project
healthy_peakset<-getPeakSet(ArchRProj = heme_archr)


for(p in unique(scATAC_ref$patient)){
  message(p)
  AML_ArchRProjects[[p]]<-addPeakSet(
    ArchRProj = AML_ArchRProjects[[p]],
    peakSet = healthy_peakset,
    genomeAnnotation = getGenomeAnnotation(AML_ArchRProjects[[p]]),
    force = T
  )
  
  AML_ArchRProjects[[p]] <- addPeakMatrix(AML_ArchRProjects[[p]])
}


## Granja et al code for remapping of disease cells onto UMAP

#Clustering and scATAC-seq UMAP for Hematopoiesis data
#06/02/19
#Cite Granja*, Klemm*, Mcginnis* et al. 
#A single cell framework for multi-omic analysis of disease identifies 
#malignant regulatory signatures in mixed phenotype acute leukemia (2019)
#Created by Jeffrey Granja
library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(uwot)
library(edgeR)
library(FNN)
library(matrixStats)
library(Rcpp)
set.seed(1)

####################################################
#Functions
####################################################

#Binarize Sparse Matrix
binarizeMat <- function(mat){
  mat@x[mat@x > 0] <- 1
  mat
}

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
  
  set.seed(1)
  
  #TF IDF LSI adapted from flyATAC
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    idx <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
    mat <- mat[idx,] 
  }else{
    idx <- which(Matrix::rowSums(mat) > 0)
    mat <- mat[idx,]
  }
  
  #Calc RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  
  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/colSm)
  idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  #Return Object
  out <- list(
    matSVD = matSVD, 
    rowSm = rowSm, 
    colSm = colSm, 
    idx = idx, 
    svd = svd, 
    binarize = binarize, 
    nComponents = nComponents,
    date = Sys.Date(),
    seed = 1)
  
  out
  
}

projectLSI <- function(mat, lsi){   
  
  #Get Same Features
  mat <- mat[lsi$idx,]
  if(lsi$binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1       
  }
  
  #Calc TF IDF
  rowsToZero <- which(lsi$rowSm == 0)
  setToZero <- which((mat@i + 1) %in% rowsToZero)
  if(length(setToZero) > 0){
    mat@x[setToZero] <- 0
  }
  
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + length(lsi$colSm) / lsi$rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  if(length(Matrix::which(is.na(tfidf),arr.ind=TRUE)) > 0){
    tfidf[Matrix::which(is.na(tfidf),arr.ind=TRUE)] <- 0 #weird Inf * 0
  }
  
  #Calc V
  V <- t(tfidf) %*% lsi$svd$u %*% diag(1/lsi$svd$d)
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svdDiag <- matrix(0, nrow=lsi$nComponents, ncol=lsi$nComponents)
  diag(svdDiag) <- lsi$svd$d
  matSVD <- t(svdDiag %*% t(V))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  return(matSVD)
  
}

#Sparse Variances Rcpp
sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)

#Compute Fast Sparse Row Variances
sparseRowVariances <- function (m){
  rM <- Matrix::rowMeans(m)
  rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
  return(rV)
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

#Clustering function using seurat SNN (Seurat v2.3.4)
seuratSNN <- function(matSVD, dims.use = 1:50, ...){
  set.seed(1)
  message("Making Seurat Object...")
  mat <- matrix(rnorm(nrow(matSVD) * 3, 1000), ncol = nrow(matSVD), nrow = 3)
  colnames(mat) <- rownames(matSVD)
  rownames(mat) <- c("1","2","3")
  obj <- Seurat::CreateSeuratObject(mat, project='scATAC', min.cells=0)
  dim_reduc<-CreateDimReducObject(embeddings=matSVD, key="pca", assay = DefaultAssay(obj))
  obj[["pca"]]<-dim_reduc
  obj <- Seurat::FindNeighbors(object = obj)
  obj <- Seurat::FindClusters(object = obj)
  # obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "cell.embeddings", new.data = matSVD)
  # obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "key", new.data = "PC")
  # obj <- Seurat::FindClusters(object = obj, reduction = "pca", dims = dims.use, verbose = TRUE, graph.name="RNA")
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  paste0("Cluster",match(clust, unique(clust)))
}


sparseMatTTest <- function(mat1, mat2, m0 = 0){
  #Get Population Values
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  n <- n1 + n2
  #Sparse Row Means
  m1 <- Matrix::rowMeans(mat1, na.rm=TRUE)
  m2 <- Matrix::rowMeans(mat2, na.rm=TRUE)
  #Sparse Row Variances
  v1 <- ArchRx:::computeSparseRowVariances(mat1@i + 1, mat1@x, m1, n1)
  v2 <- ArchRx:::computeSparseRowVariances(mat2@i + 1, mat2@x, m2, n2)
  #Calculate T Statistic
  se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2) )
  tstat <- (m1-m2-m0)/se
  #tstat <- sqrt((n1 * n2) / n) / sqrt((n1-1)/(n-2)*v1 + (n2-1)/(n-2)*v2)
  pvalue <- 2*pt(-abs(tstat), n - 2)
  fdr <- p.adjust(pvalue, method = "fdr")
  out <- data.frame(fdr = fdr, pval = pvalue, tstat = tstat, mean1 = m1, mean2 = m2, var1 = v1, var2 = v2, n1 = n1, n2 = n2)
  return(out)
}

####################################################
#Input Data
####################################################
# #Read in Summarized Experiment
# #Please Note Code here has been modified to work with finalized summarized experiment
# 
# #Reference Summarized Experiment
# #Contains Peaks for Reference Hematopoiesis only
# #seReference <- readRDS("data/Supplementary_Data_Hematopoiesis/scATAC-Healthy-Hematopoiesis-190429.rds")
# #seReference <- seReference[,sample(1:ncol(seReference),5000)] subset data to test since its faster
# 
# # substitute matrix for peaks matrix from ArchR. Subset for only healthy heme cells
# seReference<-getMatrixFromProject(
#   ArchRProj = heme_archr,
#   useMatrix = "PeakMatrix"
# )
# # dir.create("healthy_hematopoiesis_ArchR/peak_matrices/")
# # saveRDS(seReference, "healthy_hematopoiesis_ArchR/peak_matrices/normal_hematopoiesis_reference_peaks_matrix.rds")
# 
# seDisease_list<-lapply(unique(scATAC_ref$patient), function(x){
#   getMatrixFromProject(
#     ArchRProj = AML_ArchRProjects[[x]],
#     useMatrix = "PeakMatrix"
#   )
# })
# names(seDisease_list)<-unique(scATAC_ref$patient)
# 
# # for(p in unique(scATAC_ref$patient)){
# #   message(p)
# #   dir.create(paste0("AML_ArchR_Projects_hg19/",p,"/peak_matrices/"))
# #   saveRDS(seReference, paste0("AML_ArchR_Projects_hg19/",p,"/peak_matrices/normal_hematopoiesis_reference_peaks_matrix.rds"))
# # }
# 
# # run clustering and classification for all samples
# for(p in unique(scATAC_ref$patient)){
#   message(p)
#   
#   #SE Disease Cells
#   id <- p
#   #seDisease <- readRDS("../analysis/2019/Re-Analysis/Projections/ATAC/scATAC/output/Disease/MPAL1_R1/MPAL1_R1.se.rds")
#   seDisease<-seDisease_list[[p]]
#   rownames(seDisease) <- paste0(seqnames(seDisease),"_",start(seDisease),"_",end(seDisease))
#   
#   #Set Clustering Parameters
#   nPCs1 <- 1:25
#   nPCs2 <- 1:25
#   resolution <- 0.8 #clustering resolution
#   nTop <- 25000 #number of variable peaks
#   
#   #Create Matrix
#   mat <- cbind(assay(seReference), assay(seDisease))
#   
#   #Run LSI 1st Iteration
#   lsi1 <- calcLSI(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL)
#   clust1 <- seuratSNN(lsi1[[1]], dims.use = nPCs1, resolution = resolution)
#   
#   #Make Pseudo Bulk Library
#   message("Making PseudoBulk...")
#   mat <- mat[,rownames(lsi1[[1]]), drop = FALSE] #sometimes cells are filtered
#   mat@x[mat@x > 0] <- 1 #binarize
#   clusterSums <- groupSums(mat = mat, groups = clust1, sparse = TRUE) #Group Sums
#   logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
#   varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), nTop) #Top variable peaks
#   
#   #Run LSI 2nd Iteration
#   lsi2 <- calcLSI(mat[varPeaks,,drop=FALSE], nComponents = 50, binarize = TRUE, nFeatures = NULL)
#   clust2 <- seuratSNN(lsi2[[1]], dims.use = nPCs2, resolution = resolution)
#   
#   #UMAP
#   set.seed(1)
#   umap <- uwot::umap(
#     lsi2$matSVD[,1:25], 
#     n_neighbors = 55, 
#     min_dist = 0.45, 
#     metric = "euclidean", 
#     n_threads = 5, 
#     verbose = TRUE, 
#     ret_model = FALSE
#   )
#   
#   #Plot Info
#   cells <- c(
#     rep("reference", sum(rownames(lsi2$matSVD) %in% colnames(seReference))),
#     rep("disease",sum(rownames(lsi2$matSVD) %in% colnames(seDisease)))
#   )
#   
#   splitCells <- split(cells,clust2)
#   df <- data.frame(
#     clusters = names(splitCells),
#     proportion = unlist(lapply(seq_along(splitCells), function(x) sum(splitCells[[x]]=="disease") / length(splitCells[[x]])))
#   )
#   
#   #Plot UMAP Data Frame
#   plotDF <- data.frame(umap)
#   rownames(plotDF) <- c(colnames(seReference), colnames(seDisease))
#   plotDF$type <- cells
#   plotDF$clusters <- clust2
#   plotDF$classification <- 0
#   #If disease cells are clustered with healthy cluster (proportion > 0.9) we will classify these as healthy-like
#   plotDF$classification[plotDF$type == "disease" & plotDF$clusters %in% paste0(df$clusters[df[,2] > 0.9])] <- 1
#   plotDF$classification[plotDF$type == "disease"] <- plotDF$classification[plotDF$type == "disease"] + 1
#   plotDF <- plotDF[order(plotDF$classification), ]
#   
#   #Formal Classification
#   plotDF$classificationSTR <- "reference"
#   plotDF$classificationSTR[plotDF$classification==1] <- "healthy-like"
#   plotDF$classificationSTR[plotDF$classification==2] <- "disease-like"
#   
#   cellColData<-rbind(as.data.frame(heme_archr@cellColData)[,c("Sample"),drop=F],
#                      as.data.frame(AML_ArchRProjects[[p]]@cellColData)[,c("Sample"),drop=F])
#   plotDF$sample<-cellColData$Sample[match(rownames(plotDF), rownames(cellColData))]
#   
#   plotDF$sample_simple<-plotDF$sample
#   plotDF$sample_simple[plotDF$sample_simple%ni%scATAC_ref$ID[scATAC_ref$patient==p]]<-"reference"
#   
#   heme_labels<-readRDS("external_input_files/scATAC-Healthy-Hematopoiesis-191120.rds")
#   cell_labels<-data.frame(ID=paste(heme_labels$Group, heme_labels$Barcode, sep="#"),
#                           CellType=heme_labels$BioClassification)
#   plotDF$CellType<-cell_labels$CellType[match(rownames(plotDF), cell_labels$ID)]
#   
#   #Plot PDFs
#   plotDir <- paste0("outputs/scATAC_plots/scATAC_projection/")
#   dir.create(plotDir,recursive=TRUE)
#   # pdf(paste0(plotDir,id,"-Classification-UMAP.pdf"), width = 8, height = 7, useDingbats = FALSE)
#   # ggplot(plotDF, aes(X1,X2,color=classificationSTR)) +
#   #   geom_point(size=0.5) +
#   #   theme_classic() +
#   #   xlab("UMAP Dimension 1") +
#   #   ylab("UMAP Dimension 2") +
#   #   scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3"))
#   # dev.off()
#   
#   p1<-ggPoint(plotDF$X1, plotDF$X2, color = plotDF$classificationSTR, 
#               pal = c("lightgrey","dodgerblue3","firebrick3"), 
#               colorOrder = c("reference","healthy-like","disease-like"),
#               labelMeans = T, 
#               xlabel = "UMAP Dimension 1",
#               ylabel = "UMAP Dimension 2", 
#               size=0.05,
#               baseSize = 6, 
#               labelSize = 1.5)
#   
#   pdf(paste0(plotDir,id,"-Classification-UMAP.pdf"), width = 4, height = 4, useDingbats = FALSE)
#   print(p1)
#   dev.off()
#   
#   p1<-ggPoint(plotDF$X1, plotDF$X2, color = plotDF$sample_simple, 
#               pal = c("lightgrey","#272E6A","#D51F26"), 
#               colorOrder = c("reference",sort(scATAC_ref$ID[scATAC_ref$patient==p])),
#               labelMeans = T, 
#               xlabel = "UMAP Dimension 1",
#               ylabel = "UMAP Dimension 2", 
#               size=0.05,
#               baseSize = 6, 
#               labelSize = 1.5)
#   pdf(paste0(plotDir,id,"-Sample-UMAP.pdf"), width = 4, height = 4, useDingbats = FALSE)
#   print(p1)
#   dev.off()
#   
#   df_sub<-plotDF[plotDF$sample%ni%scATAC_ref$ID[scATAC_ref$patient==p],]
#   df_sub<-df_sub[!is.na(df_sub$CellType),]
#   p1<-ggPoint(df_sub$X1, df_sub$X2, color = df_sub$CellType, 
#               pal = ArchR::paletteDiscrete(unique(plotDF$CellType)), 
#               labelMeans = T, 
#               xlabel = "UMAP Dimension 1",
#               ylabel = "UMAP Dimension 2", 
#               size=0.05,
#               baseSize = 6, 
#               labelSize = 1.5)
#   
#   pdf(paste0(plotDir,id,"-reference-UMAP.pdf"), width = 6, height = 6, useDingbats = FALSE)
#   print(p1)
#   dev.off()
#   
# }


####################################################
#Project Into LSI UMAP
####################################################

#Previous Reference Summarized Experiment
#Contains Peaks for Reference Hematopoiesis only
# se <- readRDS("data/Supplementary_Data_Hematopoiesis/scATAC-Healthy-Hematopoiesis-190429.rds")
se<-readRDS("external_input_files/ArchR_Projects/healthy_hematopoiesis_plus_AML/umap_outputs/scATAC-Healthy-Hematopoiesis.rds")

#Load Saved UMAP Manifold
#umapManifold <- uwot::load_uwot("/Users/armonazizi/Bioinformatics/AML_relapse_project/scATAC/healthy_hematopoiesis_ArchR/Embeddings/Save-Uwot-UMAP-Params-IterativeLSI-fad063105669-Date-2021-06-29_Time-09-51-49.tar")
umapManifold <- uwot::load_uwot("external_input_files/ArchR_Projects/healthy_hematopoiesis_plus_AML/umap_outputs/scATAC-Hematopoiesis-UMAP-model.uwot")

# regroup cell labels and plot
refDF <- data.frame(row.names = getCellNames(heme_archr), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference")
heme_labels<-readRDS("external_input_files/scATAC-Healthy-Hematopoiesis-191120.rds")
cell_labels<-data.frame(ID=paste(heme_labels$Group, heme_labels$Barcode, sep="#"),
                        CellType=heme_labels$BioClassification,
                        CellType_original=heme_labels$BioClassification)
cell_labels$CellType[cell_labels$CellType%in%c("20_CD4.N1","22_CD4.M","24_CD8.CM","23_CD8.EM","21_CD4.N2","25_NK","19_CD8.N","13_Unk","26_Unk")]<-"T/NK-like"
cell_labels$CellType[cell_labels$CellType%in%c("17_B","16_Pre.B","18_Plasma")]<-"B/Plasma-like"
cell_labels$CellType[cell_labels$CellType%in%c("06_CLP.1","15_CLP.2")]<-"CLP-like"
cell_labels$CellType[cell_labels$CellType%in%c("07_GMP","08_GMP.Neut")]<-"GMP-like"
cell_labels$CellType[cell_labels$CellType%in%c("12_CD14.Mono.2", "11_CD14.Mono.1","10_cDC","14_Unk")]<-"Mono-like"
cell_labels$CellType[cell_labels$CellType%in%c("09_pDC","04_Early.Baso")]<-"DC/Baso-like"
cell_labels$CellType[cell_labels$CellType%in%c("02_Early.Eryth","03_Late.Eryth")]<-"Erythroid-like"
cell_labels$CellType[cell_labels$CellType%in%c("01_HSC")]<-"HSC-like"
cell_labels$CellType[cell_labels$CellType%in%c("05_CMP.LMPP")]<-"CMP/LMPP-like"
celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono-like","DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
celltype_pal_2<-ArchR::paletteDiscrete(unique(cell_labels$CellType_original))

refDF$CellType<-cell_labels$CellType[match(rownames(refDF), cell_labels$ID)]
refDF<-refDF[refDF$CellType!="Unk",]
plt<-ggPoint(refDF$X1, refDF$X2, color = refDF$CellType, 
           pal = celltype_pal, 
           labelMeans = T, 
           xlabel = "UMAP Dimension 1",
           ylabel = "UMAP Dimension 2", 
           size=0.05,
           baseSize = 6, 
           labelSize = 1.5,
           rastr = T, dpi=250)
pdf("outputs/scATAC_plots/scATAC_projection/healthy_hematopoiesis_reference_umap_reclassified_celltypes.pdf", width = 4, height = 4, useDingbats = FALSE)
print(plt)
dev.off()

# need to add peaks again
healthy_peakset<-rowRanges(se)

for(p in unique(scATAC_ref$patient)){
  message(p)
  AML_ArchRProjects[[p]]<-addPeakSet(
    ArchRProj = AML_ArchRProjects[[p]],
    peakSet = healthy_peakset,
    genomeAnnotation = getGenomeAnnotation(AML_ArchRProjects[[p]]),
    force = T
  )

  AML_ArchRProjects[[p]] <- addPeakMatrix(AML_ArchRProjects[[p]])
}

seDisease_list<-lapply(unique(scATAC_ref$patient), function(x){
  getMatrixFromProject(
    ArchRProj = AML_ArchRProjects[[x]],
    useMatrix = "PeakMatrix"
  )
})
names(seDisease_list)<-unique(scATAC_ref$patient)

# run projection and classification for all samples
for(p in unique(scATAC_ref$patient)){
  message(p)
  #SE Disease Cells
  id <- p
  #seDisease <- readRDS("../analysis/2019/Re-Analysis/Projections/ATAC/scATAC/output/Disease/MPAL1_R1/MPAL1_R1.se.rds")
  seDisease<-seDisease_list[[p]]
  rownames(seDisease) <- paste0(seqnames(seDisease),"_",start(seDisease),"_",end(seDisease))
  
  #LSI Projection Matrix
  lsiPeaks <- metadata(se)$LSIPeaks
  #lsiPeaks<-heme_archr@reducedDims$IterativeLSI_peaks$LSIFeatures
  #lsiPeaks <- paste0(lsiPeaks$seqnames,"_",lsiPeaks$start,"_",lsiPeaks$end)
  matProjectLSI <- assay(seDisease[lsiPeaks,])
  
  #LSI Project
  lsiReference <- metadata(se)$LSI
  # lsiReference <- heme_archr@reducedDims$IterativeLSI_peaks
  # lsiReference$nComponents<-lsiReference$nDimensions
  # temp_mat<-assay(seReference)
  # temp_mat@x[temp_mat@x > 0] <- 1
  # lsiReference$colSm<-Matrix::colSums(temp_mat)
  lsiProjection <- projectLSI(matProjectLSI, lsiReference)
  
  #UMAP Projection
  #Set Seed Prior to umap_transform (see uwot github)
  set.seed(1)
  umapProjection <- uwot::umap_transform(as.matrix(lsiProjection)[,1:50], umapManifold, verbose = TRUE)
  
  #Plot Projection
  refDF <- data.frame(row.names = getCellNames(heme_archr), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2])
  proDF <- data.frame(row.names = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2])
  projectionDF <- rbind(refDF, proDF)
  
  cellColData<-rbind(as.data.frame(heme_archr@cellColData)[,c("Sample"),drop=F],
                     as.data.frame(AML_ArchRProjects[[p]]@cellColData)[,c("Sample"),drop=F])
  projectionDF$Sample<-cellColData$Sample[match(rownames(projectionDF), rownames(cellColData))]
  projectionDF$Sample[projectionDF$Sample%ni%scATAC_ref$ID[scATAC_ref$patient==p]]<-"reference"
  
  plotDir <- paste0("outputs/scATAC_plots/scATAC_projection/")
  
  # pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 8, height = 7, useDingbats = FALSE)
  # ggplot(projectionDF, aes(X1,X2,color=Type)) +
  #   geom_point(size=0.5) +
  #   theme_classic() +
  #   xlab("UMAP Dimension 1") +
  #   ylab("UMAP Dimension 2") +
  #   scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3"))
  # dev.off()
  
  # p1<-ggPoint(projectionDF$X1, projectionDF$X2, color = projectionDF$Type,
  #            pal = c("lightgrey","dodgerblue3","firebrick3"),
  #            colorOrder = c("reference","healthy-like","disease-like"),
  #            labelMeans = F,
  #            xlabel = "UMAP Dimension 1",
  #            ylabel = "UMAP Dimension 2",
  #            size=0.05,
  #            baseSize = 6,
  #            labelSize = 1.5)
  # plotPDF(p1, name = paste0(id,"-Projection-UMAP.pdf"), ArchRProj = heme_archr, addDOC = FALSE, width = 4, height = 4)
  
  # pdf(paste0(plotDir,id,"-Projected-Sample-UMAP.pdf"), width = 8, height = 7, useDingbats = FALSE)
  # ggplot(projectionDF, aes(X1,X2,color=Sample)) + 
  #   geom_point(size=0.5) +
  #   theme_classic() +
  #   xlab("UMAP Dimension 1") + 
  #   ylab("UMAP Dimension 2") +
  #   scale_color_manual(values=c("reference"="lightgrey","SU892"="#D51F26","SU892B"="#272E6A"))
  # dev.off()
  
  p1<-ggPoint(projectionDF$X1, projectionDF$X2, color = projectionDF$Sample, 
              pal = c("lightgrey","#272E6A","#D51F26"), 
              colorOrder = c("reference",sort(scATAC_ref$ID[scATAC_ref$patient==p])),
              labelMeans = F, 
              xlabel = "UMAP Dimension 1",
              ylabel = "UMAP Dimension 2", 
              size=0.05,
              baseSize = 6, 
              labelSize = 1.5,
              rastr = T, dpi=250)
  pdf(paste0(plotDir,id,"-Projected-Sample-UMAP.pdf"), width = 4, height = 4, useDingbats = FALSE)
  print(p1)
  dev.off()
  
  # cell classification by nearest neighbor
  # need to do in LSI SVD space
  svdReference <- as.data.frame(lsiReference$matSVD)
  svdDisease <- as.data.frame(as.matrix(lsiProjection))
  
  #KNN Nearest Neighbor using FNN
  svdReference$CellType<-cell_labels$CellType[match(rownames(svdReference), cell_labels$ID)]
  svdReference<-svdReference[complete.cases(svdReference),]
  classifications<-class::knn(svdReference[,colnames(svdReference)!="CellType"],
                              svdDisease,
                              svdReference$CellType,
                              k = 10)
  
  #KNN Nearest Neighbor using FNN - original celltype labels
  svdReference$CellType_original<-cell_labels$CellType_original[match(rownames(svdReference), cell_labels$ID)]
  svdReference<-svdReference[complete.cases(svdReference),]
  classifications_original_labels<-class::knn(svdReference[,colnames(svdReference)%ni%c("CellType_original","CellType")],
                              svdDisease,
                              svdReference$CellType_original,
                              k = 10)
  
  
  
  # projectionDF$CellType<-cell_labels$CellType[match(rownames(projectionDF), cell_labels$ID)]
  # projectionDF_sub<-rbind(projectionDF[complete.cases(projectionDF),], projectionDF[projectionDF$Sample!="reference",])
  # classifications<-class::knn(projectionDF_sub[projectionDF_sub$Sample=="reference", c("X1","X2")],
  #                             projectionDF_sub[projectionDF_sub$Sample!="reference", c("X1","X2")],
  #                             projectionDF_sub$CellType[projectionDF_sub$Sample=="reference"],
  #                             k = 10)
  
  projectionDF$AML_classification<-as.character(classifications[match(rownames(projectionDF), rownames(svdDisease))])
  projectionDF$AML_classification[is.na(projectionDF$AML_classification)]<-"reference"
  
  projectionDF$AML_classification_original_labels<-as.character(classifications_original_labels[match(rownames(projectionDF), rownames(svdDisease))])
  projectionDF$AML_classification_original_labels[is.na(projectionDF$AML_classification_original_labels)]<-"reference"
  
  p1<-ggPoint(projectionDF$X1, projectionDF$X2, color = projectionDF$AML_classification, 
              pal = c(celltype_pal[names(celltype_pal)%in%unique(projectionDF$AML_classification)],"lightgrey"), 
              colorOrder = c(names(celltype_pal[names(celltype_pal)%in%unique(projectionDF$AML_classification)]),"reference"),
              labelMeans = F, 
              xlabel = "UMAP Dimension 1",
              ylabel = "UMAP Dimension 2", 
              size=0.05,
              baseSize = 6, 
              labelSize = 1.5,
              rastr = T, dpi=250)
  pdf(paste0(plotDir,id,"-Projected-Classification-UMAP.pdf"), width = 4, height = 4, useDingbats = FALSE)
  print(p1)
  dev.off()
  
  p1<-ggPoint(projectionDF$X1, projectionDF$X2, color = projectionDF$AML_classification_original_labels, 
              pal = c(celltype_pal_2[names(celltype_pal_2)%in%unique(projectionDF$AML_classification_original_labels)],"lightgrey"), 
              colorOrder = c(names(celltype_pal_2[names(celltype_pal_2)%in%unique(projectionDF$AML_classification_original_labels)]),"reference"),
              labelMeans = F, 
              xlabel = "UMAP Dimension 1",
              ylabel = "UMAP Dimension 2", 
              size=0.05,
              baseSize = 6, 
              labelSize = 1.5,
              rastr = T, dpi=250)
  pdf(paste0(plotDir,id,"-Projected-Classification-original-labels-UMAP.pdf"), width = 4, height = 4, useDingbats = FALSE)
  print(p1)
  dev.off()
  
  dir.create("outputs/scATAC/AML_scATAC_cell_classifications/")
  saveRDS(projectionDF, paste0("outputs/scATAC/AML_scATAC_projection/",id,"_LSI_projection.rds"))
}


################################
# Plot cell classification fractions
################################

celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono-like","DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))

for(p in unique(scATAC_ref$patient)){
  classifications<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  plot_input<-as.data.frame(table(classifications[classifications$Sample%in%scATAC_ref$ID[scATAC_ref$patient==p],c("Sample","AML_classification")]))
  for(s in scATAC_ref$ID[scATAC_ref$patient==p]){
    plot_input$Freq[plot_input$Sample==s]<-plot_input$Freq[plot_input$Sample==s]/sum(plot_input$Freq[plot_input$Sample==s])
  }
  pdf(paste0("outputs/scATAC_plots/scATAC_projection/",p,"-Classifications-Barplot.pdf"), width = 3.5, height = 4, useDingbats = FALSE)
  print(ggplot(plot_input, aes(x=Sample, y=Freq, fill=AML_classification)) +
          geom_bar(position="fill", stat="identity", color="black", size=0.2) +
          theme_classic() +
          scale_fill_manual(values = celltype_pal) +
          ggtitle(p) +
          ylab("Frequency Of Classified Cell Type") +
          theme(legend.title = element_text(size=10)))
  dev.off()
}

