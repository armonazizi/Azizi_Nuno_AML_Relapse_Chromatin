# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# 
# generate vaf plots for stanford samples


library(openxlsx)
library(dplyr)
library(tibble)
library(ComplexHeatmap)

# set to base directory
setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")

# read and format processed genotyping data
mutation_data<-openxlsx::read.xlsx("inputs/bulk_sample_info/DxRel_Variants_20200914.xlsx")
mutation_data<-mutation_data[,c("SAMPLE",  "Tier", "VarScan", "VarDict", "Mutect","Gene.refGene","Chr","Start","End","Ref","Alt", "ExonicFunc.refGene","Func.refGene")]
mutation_data$pos<-paste(mutation_data$Chr,mutation_data$Start,mutation_data$End,mutation_data$Ref,mutation_data$Alt,sep="_")
mutation_data$Case<-sapply(strsplit(mutation_data$SAMPLE,"_"), '[', 1) %>% substr(.,0,5)
mutation_data$timepoint<-sapply(strsplit(mutation_data$SAMPLE,"_"), '[', 2)
mutation_data$VarScan[mutation_data$VarScan=="."]<-NA
mutation_data$VarScan<-as.numeric(mutation_data$VarScan)
mutation_data$VAF<-rowMeans(mutation_data[,c("VarScan", "VarDict", "Mutect")], na.rm = T)


# filter
mutation_data<-mutation_data[mutation_data$ExonicFunc.refGene!="synonymous SNV",]
mutation_data<-mutation_data[mutation_data$Func.refGene!="intronic",]


pts<-c("SU142", "SU360", "SU484", "SU892")

for(p in pts){
  message(p)
  
  data<-mutation_data[mutation_data$Case==p,]
  
  combined_data<-data.frame(gene=c(),pos=c(),dx_vaf=c(),rel_vaf=c())
  
  for(pos in unique(data$pos)){
    
    if(length(data$VAF[data$pos==pos&data$timepoint=="Diagnosis"])>0){
      dx_vaf<-data$VAF[data$pos==pos&data$timepoint=="Diagnosis"]
    }else{
      dx_vaf<-0
    }
    
    if(length(data$VAF[data$pos==pos&data$timepoint=="Relapse"])>0){
      rel_vaf<-data$VAF[data$pos==pos&data$timepoint=="Relapse"]
    }else{
      rel_vaf<-0
    }
    
    group<-"Stable"
    if(dx_vaf<0.05&rel_vaf>0.1){group<-"Gain"}
    if(rel_vaf<0.05&dx_vaf>0.1){group<-"Loss"}
    
    
    combined_data<-rbind(combined_data,
                         data.frame(gene=unique(data$Gene.refGene[data$pos==pos]),
                                    pos=c(pos),
                                    dx_vaf=dx_vaf,
                                    rel_vaf=rel_vaf,
                                    group=group))
  }
  
  col = c(Stable = "#EE3A2D", Gain = "#2270B6", Loss="#70ACD4")
  
  print(combined_data)
  
  plt_input<-reshape2::melt(combined_data, id=c("gene","pos","group"))
  
  plt_input$label<-plt_input$gene
  plt_input$label[plt_input$variable=="dx_vaf"]<-""
  
  plt_input$variable<-as.character(plt_input$variable)
  
  plt_input$variable[plt_input$variable=="dx_vaf"]<-"DX VAF"
  plt_input$variable[plt_input$variable=="rel_vaf"]<-"REL VAF"
  
  plot<-ggplot(data = plt_input, aes(x = variable, y = as.numeric(value), group = pos, color=group, label=label))+
    geom_line(size=0.8) +
    geom_point() +
    theme_classic() +
    scale_color_manual(NULL,values=col) +
    scale_y_continuous(limits=c(0,1)) +
    xlab(NULL) +
    ylab("Variant Allele Frequency") +
    geom_label_repel(color="black", max.overlaps = Inf) +
    ggtitle(paste(p,"VAF Plot"))
  
  print(plot)
  
  pdf(paste0("outputs/bulk_genotyping_plots/",p,"_mutation_vaf_plot.pdf"), width = 5, height = 4)
  print(plot)
  dev.off()
  
  plot<-ggplot(data = plt_input, aes(x = variable, y = as.numeric(value), group = pos, color=group, label=label))+
    geom_line(size=0.8) +
    geom_point() +
    theme_classic() +
    scale_color_manual(NULL,values=col) +
    scale_y_continuous(limits=c(0,1)) +
    xlab(NULL) +
    ylab("Variant Allele Frequency") +
    ggtitle(paste(p,"VAF Plot"))
  
  print(plot)
  
  pdf(paste0("outputs/bulk_genotyping_plots/",p,"_mutation_vaf_plot_no_labels.pdf"), width = 5, height = 4)
  print(plot)
  dev.off()
}
