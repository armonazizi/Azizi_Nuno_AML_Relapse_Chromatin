# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# 
# analyze chromosome wide changes in bulk AML relapse samples

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")


# look at peak fold changes

diff_peaks_data<-readRDS("outputs/bulk_analysis/REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")
diff_peaks_data$chr<-sapply(strsplit(rownames(diff_peaks_data),"_"),'[',1)
diff_peaks_data$start<-as.numeric(sapply(strsplit(rownames(diff_peaks_data),"_"),'[',2))
diff_peaks_data$end<-as.numeric(sapply(strsplit(rownames(diff_peaks_data),"_"),'[',3))

for(chr in unique(diff_peaks_data$chr)){
  plt_input<-diff_peaks_data[diff_peaks_data$chr==chr,]
  
  splits<-seq(min(plt_input$start),max(plt_input$start),1000000)
  avg_data<-data.frame(start=splits,
                       avg_fc=unlist(lapply(splits, function(x) {
                         mean(plt_input$log2FoldChange[plt_input$start>(x-1000000)&plt_input$start<(x+1000000)])
                       })))
  avg_data[is.na(avg_data)]<-0
  
  
  plt<-ggplot(plt_input, aes(x=start,y=log2FoldChange)) +
    geom_point(size=0.5) +
    geom_line(data = avg_data, mapping=aes(x=start,y=avg_fc), color="red") +
    ggtitle(chr)
  
  print(plt)
}

# Analysis of genomic bin data

# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")
bulk_samples<-sample_info$ID[sample_info$is_bulk==1]
bulk_sample_info<-sample_info[sample_info$ID%in%bulk_samples,]
genotype_groups<-genotyping_info[,c("ID","genomic_bin")]
bulk_sample_info$genotype_groups<-genotype_groups$genomic_bin[match(bulk_sample_info$patient,genotype_groups$ID)]

# read ATAC count matrix
count_matrix<-readRDS("inputs/100_kb_bins_count_matrix.rds")
rownames(count_matrix)<-paste(count_matrix$Chr,count_matrix$Start,count_matrix$End,sep="_")
count_matrix<-count_matrix[,4:length(colnames(count_matrix))]

## paired DEseq analysis of stable samples
bulk_sample_info_sub<-bulk_sample_info[bulk_sample_info$genotype_groups=="Stable",]
bulk_sample_info_sub<-bulk_sample_info_sub[!is.na(bulk_sample_info_sub$name),]
deseq_matrix<-count_matrix[,bulk_sample_info_sub$name]

data_reference<-bulk_sample_info_sub[,c("name","patient","timepoint")]
colnames(data_reference)<-c("sample","patient","timepoint")
dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~timepoint+patient))

# do pca
rld <- varianceStabilizingTransformation(dds)

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
deseq_Results<-as.data.frame(deseq_Results)

deseq_Results$color<-rep("black",nrow(deseq_Results))
deseq_Results$color[deseq_Results$padj<0.05&abs(deseq_Results$log2FoldChange)>0.25]<-"red"
deseq_Results$size<-as.numeric(sapply(strsplit(rownames(deseq_Results), "_"),`[`,3))-as.numeric(sapply(strsplit(rownames(deseq_Results), "_"),`[`,2))
deseq_Results<-deseq_Results[complete.cases(deseq_Results),]

plt<-ggplot(deseq_Results, aes(x=log2FoldChange,y=-log10(padj), color=color)) +
  ggrastr::rasterise(geom_point(size=0.5), dpi=500) +
  scale_color_manual(values=c("black","darkred")) +
  theme_bw() +
  xlab("Log2FC REL - DX") +
  ylab("-log10(Adjusted p-val)") +
  ggtitle("Relapse vs Diagnosis Peaks Differential Accessibility") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color="black"))

plt

deseq_Results$chr<-sapply(strsplit(rownames(deseq_Results),"_"),'[',1)
deseq_Results$start<-as.numeric(sapply(strsplit(rownames(deseq_Results),"_"),'[',2))
deseq_Results$end<-as.numeric(sapply(strsplit(rownames(deseq_Results),"_"),'[',3))


diff_gene_data<-readRDS("outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")
gene_bodies<-as.data.frame(fread("inputs/gene_bodies/hg38_gene_bodies.bed"))
up_genes<-rownames(diff_gene_data)[diff_gene_data$log2FoldChange>0.4&diff_gene_data$padj<0.05]

up_gene_bodies<-gene_bodies[gene_bodies$gene_name%in%up_genes,]

for(chr in unique(deseq_Results$chr)){
  message(chr)
  
  plt_input<-deseq_Results[deseq_Results$chr==chr,]
  
  max_len<-max(plt_input$end)
  genes_to_plot<-up_gene_bodies[up_gene_bodies$seqnames==chr,]
  genes_to_plot$y<-rep(1,nrow(genes_to_plot))
  genes_to_plot$bin<-floor((genes_to_plot$start/max_len)*nrow(plt_input)*10)
  genes_to_plot<-genes_to_plot[!is.na(genes_to_plot$gene_name),]
  
  genes_to_plot$bin_remap<-genes_to_plot$bin
  
  if(nrow(genes_to_plot)>1){
    genes_to_plot<-genes_to_plot[order(genes_to_plot$bin),]
    for(i in 2:nrow(genes_to_plot)){
      if(genes_to_plot[i,"bin"]==genes_to_plot[i-1,"bin"]){
        genes_to_plot[i,"bin_remap"]<-genes_to_plot[i-1,"bin_remap"]+1
      }
    }
  }
  
  plt_input$moving_avg<-unlist(lapply(1:nrow(plt_input), function(x) {
    mean(plt_input$log2FoldChange[x:(x+40)])
  }))
  
  chr_plt<-ggplot(plt_input, aes(x=start,y=log2FoldChange)) +
    geom_hline(yintercept = 0, linetype = "dashed", size=1) +
    ggrastr::rasterise(geom_point(size=0.5), dpi=300) +
    geom_line(aes(x=start,y=moving_avg), color="#EE3A2D", size=1) +
    #geom_point(data = genes_to_plot, mapping=aes(x=start,y=y),color="blue") +
    theme_classic() +
    theme(title = element_text(hjust = 0.5)) +
    xlab("Chromosomal Position") +
    ylab("Accessibility Fold Change Enrichment") +
    scale_x_continuous(limits = c(0,max(plt_input$end))) +
    ggtitle(paste(chr,"Relapse vs Diagnosis Accessibility"))
  
  if(nrow(genes_to_plot)>0){
    labels<-Heatmap(t(data.frame(x=rep(0,nrow(plt_input)*10))), col = c("white"), 
                    bottom_annotation = columnAnnotation(foo = anno_mark(at = genes_to_plot$bin_remap, 
                                                                         labels = genes_to_plot$gene_name, 
                                                                         side="bottom", 
                                                                         padding = 0.6,
                                                                         labels_gp = gpar(fontsize = 8))),
                    show_heatmap_legend = F, show_row_names = F, cluster_columns = F) 
  }else{
    labels<-Heatmap(t(data.frame(x=rep(0,nrow(plt_input)))), col = c("white"), 
                    show_heatmap_legend = F, show_row_names = F, cluster_columns = F)
  }
  
  grob <- grid.grabExpr(ComplexHeatmap::draw(labels, padding = unit(c(5,23,0,8), "mm"))) 
  
  
  plt<-ggarrange(chr_plt,grob, ncol=1, heights = c(6,1))
  
  print(plt)
  
  pdf(paste0("outputs/bulk_chromosome_accessibility_plots/stable_REL_vs_DX_",chr,".pdf"), width=6, height=4)
  print(plt)
  dev.off()
}
