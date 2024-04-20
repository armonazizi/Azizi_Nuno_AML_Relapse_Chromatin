# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 2
# Stanford Relapse Samples Genotyping Oncoprint
# Armon Azizi (aazizi@stanford.edu)

library(openxlsx)
library(dplyr)
library(tibble)
library(ComplexHeatmap)

# set to base directory
setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")

# read and format processed genotyping data
genotype_data<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")
genotype_data$timepoint<-rep("DX",nrow(genotype_data))
genotype_data$timepoint[grep("\\.",genotype_data$ID)]<-"REL"
genotype_data$patient<-sapply(strsplit(genotype_data$ID, "\\."),'[',1)

patients<-unique(grep("\\.",genotype_data$ID, value=T, invert = T))

# exclude cases without data
patients<-patients[!patients%in%c("SU791","SU781","SU703")]

# get genes to plot based on recurrent genes from other studies
genes<-unique(c(sapply(strsplit(strsplit(paste(gsub(" ","",genotype_data$tier_1_variants), collapse=","),",")[[1]],"-"),'[',1),
                sapply(strsplit(strsplit(paste(gsub(" ","",genotype_data$tier_2_variants), collapse=","),",")[[1]],"-"),'[',1)))
genes<-genes[genes!="NA"]
genes<-genes[!is.na(genes)]

# parse data into matrix
oncoprint_input<-data.frame()

for(p in patients){
  oncoprint_input[,as.character(p)]<-rep(NA,nrow(oncoprint_input))
  dx_genes<-unique(c(sapply(strsplit(strsplit(gsub(" ","",genotype_data[genotype_data$patient==p&genotype_data$timepoint=="DX","tier_1_variants"]),",")[[1]],"-"),'[',1),
                     sapply(strsplit(strsplit(gsub(" ","",genotype_data[genotype_data$patient==p&genotype_data$timepoint=="DX","tier_2_variants"]),",")[[1]],"-"),'[',1)))
  dx_genes<-dx_genes[dx_genes!="NA"]
  dx_genes<-dx_genes[!is.na(dx_genes)]
  
  rel_genes<-unique(c(sapply(strsplit(strsplit(gsub(" ","",genotype_data[genotype_data$patient==p&genotype_data$timepoint=="REL","tier_1_variants"]),",")[[1]],"-"),'[',1),
                     sapply(strsplit(strsplit(gsub(" ","",genotype_data[genotype_data$patient==p&genotype_data$timepoint=="REL","tier_2_variants"]),",")[[1]],"-"),'[',1)))
  rel_genes<-rel_genes[rel_genes!="NA"]
  rel_genes<-rel_genes[!is.na(rel_genes)]
  
  all_genes<-unique(c(dx_genes,rel_genes))
  
  for(g in all_genes){
    print(g)
    print(p)
    
    if((g%in%dx_genes) & (g%in%rel_genes)){
      oncoprint_input[g,as.character(p)]<-3
    }else if((g%in%dx_genes) & !(g%in%rel_genes)){
      oncoprint_input[g,as.character(p)]<-1
    }else if(!(g%in%dx_genes) & (g%in%rel_genes)){
      oncoprint_input[g,as.character(p)]<-2
    }
  }
}
oncoprint_input[is.na(oncoprint_input)]<-0

# format oncoprint input
oncoprint_input<-oncoprint_input[rev(order(rowSums(oncoprint_input>0))),]

oncoprint_input<-as.data.frame(t(oncoprint_input))
oncoprint_input<-rownames_to_column(oncoprint_input)
oncoprint_input<-arrange_at(oncoprint_input, colnames(oncoprint_input)[2:ncol(oncoprint_input)])
oncoprint_input<-column_to_rownames(oncoprint_input)
oncoprint_input<-t(oncoprint_input)
oncoprint_input<-oncoprint_input[,rev(colnames(oncoprint_input))]
oncoprint_input<-as.data.frame(oncoprint_input)
oncoprint_input[] <- lapply(oncoprint_input, as.character)
oncoprint_input[oncoprint_input==0]<-""
oncoprint_input[oncoprint_input==1]<-"Lost"
oncoprint_input[oncoprint_input==2]<-"Gained"
oncoprint_input[oncoprint_input==3]<-"Stable"

# parse FLT3, NPM1, and karyotype calls
flt3_mut<-c()
for(p in colnames(oncoprint_input)){
  dx_mut<-genotype_data$FLT3_mut[genotype_data$patient==p&genotype_data$timepoint=="DX"]
  rel_mut<-genotype_data$FLT3_mut[genotype_data$patient==p&genotype_data$timepoint=="REL"]
  
  if(dx_mut==1&rel_mut==1){
    flt3_mut<-c(flt3_mut,"Stable")
  }else if(dx_mut==1&rel_mut!=1){
    flt3_mut<-c(flt3_mut,"Lost")
  }else if(dx_mut!=1&rel_mut==1){
    flt3_mut<-c(flt3_mut,"Gained")
  }else if(dx_mut=="TBD"|rel_mut=="TBD"){
    flt3_mut<-c(flt3_mut,"TBD")
  }else{
    flt3_mut<-c(flt3_mut,"")
  }
}

npm1_mut<-c()
for(p in colnames(oncoprint_input)){
  dx_mut<-genotype_data$NPM1_mut[genotype_data$patient==p&genotype_data$timepoint=="DX"]
  rel_mut<-genotype_data$NPM1_mut[genotype_data$patient==p&genotype_data$timepoint=="REL"]
  
  if(dx_mut==1&rel_mut==1){
    npm1_mut<-c(npm1_mut,"Stable")
  }else if(dx_mut==1&rel_mut!=1){
    npm1_mut<-c(npm1_mut,"Lost")
  }else if(dx_mut!=1&rel_mut==1){
    npm1_mut<-c(npm1_mut,"Gained")
  }else if(dx_mut=="TBD"|rel_mut=="TBD"){
    npm1_mut<-c(npm1_mut,"TBD")
  }else{
    npm1_mut<-c(npm1_mut,"")
  }
}

karyotype<-c()
for(p in colnames(oncoprint_input)){
  dx_mut<-genotype_data$Karyotype[genotype_data$patient==p&genotype_data$timepoint=="DX"]
  rel_mut<-genotype_data$Karyotype[genotype_data$patient==p&genotype_data$timepoint=="REL"]
  
  if(dx_mut=="Complex"&rel_mut=="Complex"){
    karyotype<-c(karyotype,"Stable")
  }else if(dx_mut=="Complex"&rel_mut=="NK"){
    karyotype<-c(karyotype,"Lost")
  }else if(dx_mut=="NK"&rel_mut=="Complex"){
    karyotype<-c(karyotype,"Gained")
  }else if(dx_mut=="N/A"|rel_mut=="N/A"){
    karyotype<-c(karyotype,"Unknown")
  }else{
    karyotype<-c(karyotype,"")
  }
}





col = c(Stable = "#EE3A2D", Gained = "#2270B6", Lost="#70ACD4") ####changed colors

op<-oncoPrint(as.matrix(oncoprint_input),
          alter_fun = function(x, y, w, h, v){
            n = sum(v)
            h = h*0.9
            if(n){
              grid.rect(x, y - h*0.5 + 1:n/n*h, w*1, 1/n*h,gp = gpar(fill = col[names(which(v))], col = "white"), just = "top")} ###changed heights slightly and added white border
            else{grid.rect(x, y, w, h,gp = gpar(fill = "#E9E9E9", col = "white"))}}, ###adds grey background if field is empty
          col = col,
          row_order = NULL,
          column_order = colnames(oncoprint_input),
          row_names_side="left",
          #row_names_gp = gpar(fontsize=5),
          show_column_names = T,
          pct_side="right",
          #pct_gp = gpar(fontsize=5),
          remove_empty_columns=TRUE)

# save input
#saveRDS(oncoprint_input, "oncoprint.rds")
# 
# # Plotting and save as pdf
# op<-oncoPrint(as.matrix(oncoprint_input),
#           alter_fun = list(
#             background = function(x, y, w, h) grid.rect(x, y, w, h, 
#                                                         gp = gpar(fill = "grey90")),
#             Stable = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.8, 
#                                                  gp = gpar(fill = "red", col = NA)),
#             Gained = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.8, 
#                                                     gp = gpar(fill = "blue", col = NA)),
#             Lost = function(x, y, w, h) grid.rect(x, y, w*0.8, h*0.8, 
#                                                     gp = gpar(fill = "purple", col = NA))
#           ),
#           col=c(Stable="red",Gained="blue",Lost="purple"),
#           row_names_side="left",
#           row_names_gp = gpar(fontsize=5),
#           show_column_names = T,
#           pct_side="right",
#           pct_gp = gpar(fontsize=5),
#           row_order = NULL,
#           column_order = colnames(oncoprint_input))
#           # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
#           #                                    FLT3_Clinical = flt3_mut,
#           #                                    NPM1_Clinical = npm1_mut,
#           #                                    complex_karyotype = karyotype,
#           #                                    col=list(FLT3_Clinical=c(Stable="darkblue",Gained="darkred",Lost="goldenrod",TBD="black"),
#           #                                             NPM1_Clinical=c(Stable="darkblue",Gained="darkred",Lost="goldenrod",TBD="black"),
#           #                                             complex_karyotype=c(Stable="darkblue",Gained="darkred",Lost="goldenrod",Unknown="black")),
#           #                                    annotation_name_side="left"))
# 

pdf("outputs/figure_2/oncoprint.pdf", width = 6.5, height = 6)
print(op)
dev.off()
