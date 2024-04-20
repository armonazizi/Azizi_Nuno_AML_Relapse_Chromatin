# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 2
# Mapping of bulk AML ATAC samples to healthy hematopoiesis single cell manifold
#
# Much of the projection code is borrowed and modified from:
# Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. 
# Nature Biotechnology (Granja JM*, Klemm SK*, McGinnis LM*, et al. 2019)
#
# https://github.com/GreenleafLab/MPAL-Single-Cell-2019
#
# Workflow for this script is as follows:
# 1. Define functions
# 2. Project healthy celltypes to healthy single cell manifold (for validation)
# 3. Project bulk AMLs to healthy single cell manifold
# 4. Analyze projections and do nearest celltype analysis

library(ArchR)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

#### Define Functions ####
# (borrowed and modified from Granja/Klemm et al. NBiotech 2019)
safelapply <- function(..., threads = 1, preschedule = FALSE){
  
  if(tolower(.Platform$OS.type) == "windows"){
    threads <- 1
  }
  
  if(threads > 1){
    
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)
    
    errorMsg <- list()
    
    for(i in seq_along(o)){ #Make Sure this doesnt explode!
      if(inherits(o[[i]], "try-error")){
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ", i, " : "), capOut), "\n")
      }
    }
    
    if(length(errorMsg) != 0){
      
      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)
      
    }
    
  }else{
    
    o <- lapply(...)
    
  }
  
  o
  
}

safeSubset <- function(mat = NULL, subsetRows = NULL, subsetCols = NULL){
  
  if(!is.null(subsetRows)){
    idxNotIn <- which(subsetRows %ni% rownames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(length(idxNotIn), ncol = ncol(mat)))
      rownames(matNotIn) <- subsetNamesNotIn
      mat <- rbind(mat, matNotIn)
    }
    mat <- mat[subsetRows,]
  }
  
  if(!is.null(subsetCols)){
    idxNotIn <- which(subsetCols %ni% colnames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetCols[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(nrow(mat), ncol = length(idxNotIn)))
      colnames(matNotIn) <- subsetNamesNotIn
      mat <- cbind(mat, matNotIn)
    }
    mat <- mat[,subsetCols]
  }
  
  mat
  
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

# Function to simulate single cell data from bulk atac and project to a preexisting umap manifold
project_bulk<-function(lsiReference=NULL,
                       lsiPeaks=NULL,
                       seATAC=NULL,
                       umapManifold=NULL,
                       n = 250){
  
  # Reduced Dimensions
  
  rD <- lsiReference
  
  rDFeatures <- lsiPeaks
  rDGR<-rDFeatures
  
  subATAC <- subsetByOverlaps(seATAC, rDGR, ignore.strand = TRUE)
  subATAC <- subATAC[order(rowSums(as.matrix(SummarizedExperiment::assays(subATAC)[["counts"]])), decreasing = TRUE), ]
  o <- DataFrame(findOverlaps(subATAC, rDGR, ignore.strand = TRUE))
  sumOverlap <- length(unique(o[,2]))
  
  message(paste0("Overlap Ratio of Reduced Dims Features = ", (sumOverlap / length(rDGR))))
  
  o <- o[!duplicated(o$subjectHits),]
  subATAC <- subATAC[o$queryHits, ]
  rownames(subATAC) <- paste0("f", o$subjectHits)
  
  
  # Create Bulk Matrix
  
  bulkMat <- safeSubset(
    mat = SummarizedExperiment::assays(subATAC)[["counts"]], 
    subsetRows = paste0("f", seq_along(rDGR))
  )
  
  # Simulate and Project
  depthN <- round(sum(rD$rowSm / length(rD$rowSm)))
  nRep <- 5
  n2 <- ceiling(n / nRep)
  ratios <- c(2, 1.5, 1, 0.5, 0.25) #range of ratios of number of fragments
  
  simRD <- safelapply(seq_len(ncol(bulkMat)), function(x){
    counts <- bulkMat[, x]
    counts <- rep(seq_along(counts), counts)
    simMat <- lapply(seq_len(nRep), function(y){
      ratio <- ratios[y]
      simMat <- matrix(sample(x = counts, size = ceiling(ratio * depthN) * n2, replace = TRUE), ncol = n2)
      simMat <- Matrix::summary(as(simMat, "dgCMatrix"))[,-1,drop=FALSE]
      simMat[,1] <- simMat[,1] + (y - 1) * n2
      simMat
    }) %>%  Reduce("rbind", .)
    simMat <- Matrix::sparseMatrix(i = simMat[,2], j = simMat[,1], x = rep(1, nrow(simMat)), dims = c(length(rDGR), n2 * nRep))
    simRD <- as.matrix(projectLSI(simMat, rD))
    rownames(simRD) <- paste0(colnames(bulkMat)[x], "#", seq_len(nrow(simRD)))
    simRD
  }, threads = 4) %>% Reduce("rbind", .)
  
  
  umapProjection <- uwot::umap_transform(as.matrix(simRD), umapManifold, verbose = TRUE)
  rownames(umapProjection)<-rownames(simRD)
  list(simRD, umapProjection)
}


#### Analysis Start ####

# Load common objects
heme_archr<-loadArchRProject("external_input_files/ArchR_Projects/healthy_hematopoiesis_ArchR")
se<-readRDS("external_input_files/ArchR_Projects/healthy_hematopoiesis_plus_AML/umap_outputs/scATAC-Healthy-Hematopoiesis.rds")
umapManifold <- uwot::load_uwot("external_input_files/ArchR_Projects/healthy_hematopoiesis_plus_AML/umap_outputs/scATAC-Hematopoiesis-UMAP-model.uwot")
heme_labels<-readRDS("external_input_files/scATAC-Healthy-Hematopoiesis-191120.rds")


#### Healthy heme projection ####

# read bulk ATAC data. start with healthy heme.
bulk_ATAC_counts<-as.data.frame(fread("external_input_files/healthy_hematopoiesis_ATAC/ATACseq_All_Counts.txt"))

# make summarized experiment of bulk data
bulk_ATAC_counts_sub<-bulk_ATAC_counts[,grep("SU", colnames(bulk_ATAC_counts), invert = T)]
bulk_ATAC_counts_sub<-bulk_ATAC_counts_sub[,c(4:ncol(bulk_ATAC_counts_sub))]
bulk_ATAC_peaks<-GRanges(seqnames=bulk_ATAC_counts$Chr,
                         ranges = IRanges(start=bulk_ATAC_counts$Start, end=bulk_ATAC_counts$End))
bulk_ATAC_se<-SummarizedExperiment(assays = list(counts=as.matrix(bulk_ATAC_counts_sub)),
                                   rowRanges = bulk_ATAC_peaks)

# project healthy ATAC to healthy heme manifold
umapProjection<-project_bulk(lsiReference=metadata(se)$LSI,
                             lsiPeaks=rowRanges(se)[metadata(se)$LSIPeaks],
                             seATAC=bulk_ATAC_se,
                             umapManifold=umapManifold,
                             n = 250)
umapProjection<-umapProjection[[2]]

refDF <- data.frame(row.names = getCellNames(heme_archr), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2])
proDF <- data.frame(row.names = rownames(umapProjection), X1 = umapProjection[,1], X2 = umapProjection[,2])
proDF$sample<-strsplit(rownames(proDF),"#") %>% sapply('[',1)
proDF$donor<-strsplit(proDF$sample, "-") %>% sapply('[',1)

# make annotations
bulk_annotations<-as.data.frame(fread("external_input_files/healthy_hematopoiesis_ATAC/sample_annotations.txt"))
colnames(bulk_annotations)<-c("Sample", "annotation")
proDF$annotation<-bulk_annotations$annotation[match(as.vector(proDF$sample), bulk_annotations$Sample)]
refDF$sample<-rep("scATAC", nrow(as.data.frame(refDF)))
refDF$donor<-rep("scATAC", nrow(as.data.frame(refDF)))
refDF$annotation<-rep("scATAC", nrow(as.data.frame(refDF)))
projected_bulk_umap<-as.data.frame(rbind(refDF, proDF))

# label celltypes
# aggregate celltypes and replot
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
celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono/DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))

projected_bulk_umap$CellType<-cell_labels$CellType[match(rownames(projected_bulk_umap), cell_labels$ID)]

# do plotting
custom_pal<-ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))
custom_pal["scATAC"]<-"lightgrey"
p1<-ggPoint(projected_bulk_umap$X1, projected_bulk_umap$X2, color = projected_bulk_umap$annotation, 
            pal=custom_pal,
            #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
            #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
            labelMeans = T, 
            xlabel = "UMAP Dimension 1",
            ylabel = "UMAP Dimension 2", 
            size=0.05, baseSize = 8, labelSize = 2,
            rastr = T, dpi = 300)

custom_pal<-ArchR::paletteDiscrete(unique(projected_bulk_umap$donor))
custom_pal["scATAC"]<-"lightgrey"
p2<-ggPoint(projected_bulk_umap$X1, projected_bulk_umap$X2, color = projected_bulk_umap$donor, 
            pal=custom_pal,
            #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
            #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
            labelMeans = T, 
            xlabel = "UMAP Dimension 1",
            ylabel = "UMAP Dimension 2", 
            size=0.05, baseSize = 8, labelSize = 2,
            rastr = T, dpi = 300)

projected_bulk_umap_sub<-projected_bulk_umap[projected_bulk_umap$annotation=="scATAC",]
projected_bulk_umap_sub<-projected_bulk_umap_sub[complete.cases(projected_bulk_umap_sub),]
p3<-ggPoint(projected_bulk_umap_sub$X1, projected_bulk_umap_sub$X2, color = projected_bulk_umap_sub$CellType, 
            pal=custom_pal,
            #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
            #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
            labelMeans = T, 
            xlabel = "UMAP Dimension 1",
            ylabel = "UMAP Dimension 2", 
            size=0.05, baseSize = 8, labelSize = 2,
            rastr = T, dpi = 300)

plt<-ggAlignPlots(p3,p1,p2, type = "h", draw = F)

pdf("outputs/bulk_ATAC_projection_plots/healthy_hematopoiesis_bulk_projection.pdf", width = 12, height = 6)
grid::grid.draw(plt)
dev.off()


#### Bulk AML Projection ####

# bulk AML
bulk_ATAC_counts<-readRDS("external_input_files/AML_ATAC/count_matrix_raw_hg19.rds")
bulk_ATAC_peaks<-GRanges(seqnames=bulk_ATAC_counts$Chr,
                         ranges = IRanges(start=bulk_ATAC_counts$Start, end=bulk_ATAC_counts$End))
bulk_ATAC_counts<-bulk_ATAC_counts[,c(4:ncol(bulk_ATAC_counts))]
bulk_ATAC_se<-SummarizedExperiment(assays = list(counts=as.matrix(bulk_ATAC_counts)),
                                   rowRanges = bulk_ATAC_peaks)

bulkProjection<-project_bulk(lsiReference=metadata(se)$LSI,
                             lsiPeaks=rowRanges(se)[metadata(se)$LSIPeaks],
                             seATAC=bulk_ATAC_se,
                             umapManifold=umapManifold,
                             n = 250)

umapProjection<-bulkProjection[[2]]

refDF <- data.frame(row.names = getCellNames(heme_archr), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2])
proDF <- data.frame(row.names = rownames(umapProjection), X1 = umapProjection[,1], X2 = umapProjection[,2])
refDF$sample<-strsplit(rownames(refDF),"#") %>% sapply('[',1)
proDF$sample<-strsplit(rownames(proDF),"#") %>% sapply('[',1)


# make annotations
if(T){
  bulk_annotations<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
  proDF$patient<-bulk_annotations$patient[match(proDF$sample, bulk_annotations$name)]
  proDF$type<-bulk_annotations$type[match(proDF$sample, bulk_annotations$name)]
  proDF$timepoint<-bulk_annotations$timepoint[match(proDF$sample, bulk_annotations$name)]
  refDF$patient<-rep("scATAC", nrow(as.data.frame(refDF)))
  refDF$type<-rep("scATAC", nrow(as.data.frame(refDF)))
  refDF$timepoint<-rep("scATAC", nrow(as.data.frame(refDF)))
  projected_bulk_umap<-as.data.frame(rbind(refDF, proDF))
  projected_bulk_umap<-projected_bulk_umap[complete.cases(projected_bulk_umap),]
  
  # label celltypes
  # aggregate celltypes and replot
  cell_labels<-data.frame(ID=paste(heme_labels$Group, heme_labels$Barcode, sep="#"),
                          CellType=heme_labels$BioClassification)
  cell_labels$CellType[cell_labels$CellType%in%c("20_CD4.N1","22_CD4.M","24_CD8.CM","23_CD8.EM","21_CD4.N2","25_NK","19_CD8.N","13_Unk","26_Unk")]<-"T/NK-like"
  cell_labels$CellType[cell_labels$CellType%in%c("17_B","16_Pre.B","18_Plasma")]<-"B/Plasma-like"
  cell_labels$CellType[cell_labels$CellType%in%c("06_CLP.1","15_CLP.2")]<-"CLP-like"
  cell_labels$CellType[cell_labels$CellType%in%c("07_GMP","08_GMP.Neut")]<-"GMP-like"
  cell_labels$CellType[cell_labels$CellType%in%c("12_CD14.Mono.2", "11_CD14.Mono.1","10_cDC","14_Unk")]<-"Mono-like"
  cell_labels$CellType[cell_labels$CellType%in%c("09_pDC","04_Early.Baso")]<-"DC/Baso-like"
  cell_labels$CellType[cell_labels$CellType%in%c("02_Early.Eryth","03_Late.Eryth")]<-"Erythroid-like"
  cell_labels$CellType[cell_labels$CellType%in%c("01_HSC")]<-"HSC-like"
  cell_labels$CellType[cell_labels$CellType%in%c("05_CMP.LMPP")]<-"CMP/LMPP-like"
  celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono/DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
  
  projected_bulk_umap$CellType<-cell_labels$CellType[match(rownames(projected_bulk_umap), cell_labels$ID)]
}

# do plotting
custom_pal<-ArchR::paletteDiscrete(unique(projected_bulk_umap$type))
custom_pal["scATAC"]<-"lightgrey"
p1<-ggPoint(projected_bulk_umap$X1, projected_bulk_umap$X2, color = projected_bulk_umap$type, 
            pal=custom_pal,
            #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
            #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
            labelMeans = T, 
            xlabel = "UMAP Dimension 1",
            ylabel = "UMAP Dimension 2", 
            size=0.05, baseSize = 8, labelSize = 2,
            rastr = T, dpi = 300)

custom_pal<-c("#D51F26","#272E6A","lightgrey")
names(custom_pal)<-c("REL","DX","scATAC")
p2<-ggPoint(projected_bulk_umap$X1, projected_bulk_umap$X2, color = projected_bulk_umap$timepoint, 
            pal=custom_pal,
            #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
            #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
            labelMeans = T, 
            xlabel = "UMAP Dimension 1",
            ylabel = "UMAP Dimension 2", 
            size=0.05, baseSize = 8, labelSize = 2,
            rastr = T, dpi = 300)

projected_bulk_umap_sub<-projected_bulk_umap[projected_bulk_umap$type=="scATAC",]
projected_bulk_umap_sub<-projected_bulk_umap_sub[complete.cases(projected_bulk_umap_sub),]
custom_pal<-ArchR::paletteDiscrete(unique(projected_bulk_umap_sub$CellType))
custom_pal["scATAC"]<-"lightgrey"
p3<-ggPoint(projected_bulk_umap_sub$X1, projected_bulk_umap_sub$X2, color = projected_bulk_umap_sub$CellType, 
            pal=custom_pal,
            #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
            #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
            labelMeans = T, 
            xlabel = "UMAP Dimension 1",
            ylabel = "UMAP Dimension 2", 
            size=0.05, baseSize = 8, labelSize = 2,
            rastr = T, dpi = 300)


plt<-ggAlignPlots(p3,p1,p2, type = "h", draw = F)

pdf("outputs/bulk_ATAC_projection_plots/all_AML_bulk_projection.pdf", width = 12, height = 6)
grid::grid.draw(plt)
dev.off()

# separate plotting for each patient

#projected_bulk_umap$timepoint_type<-paste(projected_bulk_umap$type, projected_bulk_umap$timepoint, sep="_")
for(p in unique(projected_bulk_umap$patient)){
  if(p=="scATAC"){next}
  message(paste("Plotting", p))
  projected_bulk_umap_sub<-projected_bulk_umap[projected_bulk_umap$patient%in%c(p,"scATAC"),]
  plot_list<-list()
  for(t in unique(projected_bulk_umap_sub$type)){
    if(t!="scATAC"){
      projected_bulk_umap_sub_type<-projected_bulk_umap_sub[projected_bulk_umap_sub$type%in%c(t,"scATAC"),]
      custom_pal<-c("#D51F26","#272E6A","lightgrey")
      names(custom_pal)<-c("REL","DX","scATAC")
      p1<-ggPoint(projected_bulk_umap_sub_type$X1, projected_bulk_umap_sub_type$X2, color = projected_bulk_umap_sub_type$timepoint, 
                  pal=custom_pal,
                  labelMeans = T, 
                  xlabel = "UMAP Dimension 1",
                  ylabel = "UMAP Dimension 2", 
                  title = paste(p,t),
                  size=0.05, baseSize = 8, labelSize = 2,
                  rastr = T, dpi = 300)
      plot_list[[t]]<-p1
    }
  }
  plt<-ggAlignPlots(plotList = plot_list, type = "h", draw = F)
  
  pdf(paste0("outputs/bulk_ATAC_projection_plots/",p,"_bulk_projection_umap.pdf"), width = 4*length(plot_list), height = 4)
  grid::grid.draw(plt)
  dev.off()
}


# Plot mean projection of all bulk samples
healthy_atac_umap<-projected_bulk_umap[projected_bulk_umap$patient=="scATAC",]
healthy_atac_umap<-healthy_atac_umap[complete.cases(healthy_atac_umap),]
projected_bulk_umap_mean<-projected_bulk_umap[projected_bulk_umap$patient!="scATAC",]
projected_bulk_umap_mean<-aggregate(.~patient+type+timepoint, projected_bulk_umap_mean[,c("X1","X2","patient","timepoint","type")], mean)
projected_bulk_umap_mean<-projected_bulk_umap_mean[order(projected_bulk_umap_mean$timepoint),]

# filter to only the bulk samples
bulk_sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
bulk_sample_info<-bulk_sample_info[bulk_sample_info$is_bulk==1,]

projected_bulk_umap_mean<-projected_bulk_umap_mean[paste(projected_bulk_umap_mean$patient,projected_bulk_umap_mean$type,projected_bulk_umap_mean$timepoint,sep="_")%in%paste(bulk_sample_info$patient,bulk_sample_info$type,bulk_sample_info$timepoint,sep="_"),]


custom_pal<-ArchR::paletteDiscrete(unique(healthy_atac_umap$CellType))
custom_pal["scATAC"]<-"lightgrey"
plt<-ggPoint(healthy_atac_umap$X1, healthy_atac_umap$X2, color = healthy_atac_umap$CellType, 
             pal=custom_pal,
             #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
             #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
             labelMeans = T, 
             xlabel = "UMAP Dimension 1",
             ylabel = "UMAP Dimension 2", 
             size=0.05, baseSize = 8, labelSize = 2,
             rastr = T, dpi = 300) +
  #geom_point(projected_bulk_umap_mean[projected_bulk_umap_mean$type=="Blast",], mapping = aes(x=X1,y=X2, group=patient, color=timepoint)) +
  geom_path(projected_bulk_umap_mean, 
            mapping = aes(x=X1,y=X2, group=patient), 
            color="black", 
            arrow = arrow(length=unit(0.2,"cm"), type="open"),
            size=0.5)

pdf("outputs/bulk_ATAC_projection_plots/bulk_AML_mean_projection.pdf", width = 6, height = 6)
plt
dev.off()



# Plot mean projection of all LSC samples
healthy_atac_umap<-projected_bulk_umap[projected_bulk_umap$patient=="scATAC",]
healthy_atac_umap<-healthy_atac_umap[complete.cases(healthy_atac_umap),]
projected_bulk_umap_mean<-projected_bulk_umap[projected_bulk_umap$patient!="scATAC",]
projected_bulk_umap_mean<-aggregate(.~patient+type+timepoint, projected_bulk_umap_mean[,c("X1","X2","patient","timepoint","type")], mean)
projected_bulk_umap_mean<-projected_bulk_umap_mean[order(projected_bulk_umap_mean$timepoint),]

# filter to only the bulk samples
bulk_sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
LSC_sample_info<-bulk_sample_info[bulk_sample_info$type=="LSC",]

projected_bulk_umap_mean<-projected_bulk_umap_mean[paste(projected_bulk_umap_mean$patient,projected_bulk_umap_mean$type,projected_bulk_umap_mean$timepoint,sep="_")%in%paste(LSC_sample_info$patient,LSC_sample_info$type,LSC_sample_info$timepoint,sep="_"),]


custom_pal<-ArchR::paletteDiscrete(unique(healthy_atac_umap$CellType))
custom_pal["scATAC"]<-"lightgrey"


plt<-ggPoint(healthy_atac_umap$X1, healthy_atac_umap$X2, color = healthy_atac_umap$CellType, 
             pal=custom_pal,
             #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
             #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
             labelMeans = F, 
             xlabel = "UMAP Dimension 1",
             ylabel = "UMAP Dimension 2", 
             size=0.05, baseSize = 8, labelSize = 2,
             rastr = T, dpi = 300) +
  #geom_point(projected_bulk_umap_mean[projected_bulk_umap_mean$type=="Blast",], mapping = aes(x=X1,y=X2, group=patient, color=timepoint)) +
  geom_path(projected_bulk_umap_mean, 
            mapping = aes(x=X1,y=X2, group=patient), 
            color="black", 
            arrow = arrow(length=unit(0.2,"cm"), type="open"),
            size=0.5)

pdf("outputs/bulk_ATAC_projection_plots/LSC_AML_mean_projection.pdf", width = 6, height = 6)
plt
dev.off()




plt<-ggPoint(healthy_atac_umap$X1, healthy_atac_umap$X2, color = healthy_atac_umap$CellType, 
             pal=custom_pal,
             #pal = c(ArchR::paletteDiscrete(unique(projected_bulk_umap$annotation))), 
             #colorOrder = c("scATAC",unique(projected_bulk_umap$annotation)),
             labelMeans = F, 
             xlabel = "UMAP Dimension 1",
             ylabel = "UMAP Dimension 2", 
             size=0.5, baseSize = 8, labelSize = 2,
             rastr = T, dpi = 300, xlim = c(-5,5), ylim=c(-7,2)) +
  #geom_point(projected_bulk_umap_mean[projected_bulk_umap_mean$type=="Blast",], mapping = aes(x=X1,y=X2, group=patient, color=timepoint)) +
  geom_path(projected_bulk_umap_mean, 
            mapping = aes(x=X1,y=X2, group=patient), 
            color="black", 
            arrow = arrow(length=unit(0.3,"cm"), type="open"),
            size=0.75)

pdf("outputs/bulk_ATAC_projection_plots/LSC_AML_mean_projection_zoomed.pdf", width = 6, height = 6)
plt
dev.off()



#### Nearest celltype classification ####

## classification of bulk samples and save
# cell classification by nearest neighbor
# need to do in LSI SVD space
svdReference <- as.data.frame(metadata(se)$LSI$matSVD)
svdDisease <- as.data.frame(as.matrix(bulkProjection[[1]]))

#KNN Nearest Neighbor using FNN
svdReference$CellType<-cell_labels$CellType[match(rownames(svdReference), cell_labels$ID)]
svdReference<-svdReference[complete.cases(svdReference),]
classifications<-class::knn(svdReference[,colnames(svdReference)!="CellType"],
                            svdDisease,
                            svdReference$CellType,
                            k = 10)


projected_bulk_umap$AML_classification<-as.character(classifications[match(rownames(projected_bulk_umap), rownames(svdDisease))])
saveRDS(projected_bulk_umap, "outputs/bulk_analysis/bulk_pseudo_sc_KNN_classifications.rds")




#### Quantification of blast nearest celltype ####
projected_bulk_umap<-readRDS("outputs/bulk_analysis/bulk_pseudo_sc_KNN_classifications.rds")

# find the most common celltype for each bulk sample
library(dplyr)
most_common_celltype <- projected_bulk_umap[,c("sample", "patient","type", "timepoint", "CellType", "AML_classification")] %>%
  dplyr::group_by(sample, type, timepoint, AML_classification) %>%
  dplyr::count() %>%
  dplyr::group_by(sample, type, timepoint) %>%
  dplyr::top_n(1, n) %>%
  dplyr::ungroup()

most_common_celltype<-most_common_celltype[most_common_celltype$timepoint!="scATAC",]
most_common_celltype<-most_common_celltype[most_common_celltype$type=="Blast",]

celltype_pal<-ArchR::paletteDiscrete(c("B/Plasma-like","Mono-like", "GMP-like", "T/NK-like","DC/Baso-like","CMP/LMPP-like","CLP-like","Erythroid-like","HSC-like"))

# Create stacked barplot with p value
fisher_result <- fisher.test(table(most_common_celltype$AML_classification, most_common_celltype$timepoint))

p_value <- format(fisher_result$p.value, digits = 2)
plt<-ggplot(most_common_celltype, aes(x = timepoint, fill = AML_classification)) +
  geom_bar(position = position_fill(), color="black", linewidth=0.5) +
  scale_fill_manual(values = celltype_pal) +
  theme_classic() +
  labs(title = "Cell Type Distribution by Timepoint",
       x = "Timepoint",
       y = "Fraction Of Samples",
       fill = "AML Classification",
       subtitle = paste0("Fisher's exact test p-value: ", p_value))

print(plt)

pdf("outputs/bulk_ATAC_projection_plots/DX_vs_REL_all_bulk_AML_celltype_fractions.pdf", width = 4, height = 4)
plt
dev.off()
