# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 2
# Differential accessibility analysis between diagnosis and relapse samples for stable clonality and unstable clonality groups.

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")
bulk_samples<-sample_info$ID[sample_info$is_bulk==1]
bulk_sample_info<-sample_info[sample_info$ID%in%bulk_samples,]
genotype_groups<-genotyping_info[,c("ID","genomic_bin")]
bulk_sample_info$genotype_groups<-genotype_groups$genomic_bin[match(bulk_sample_info$patient,genotype_groups$ID)]


# plot fraction of each clonal group
plt_input<-genotyping_info[genotyping_info$ID%in%sample_info$patient,c("ID","genomic_bin")]
plt_input<-plt_input[complete.cases(plt_input),]
plt_input$genomic_bin[plt_input$genomic_bin=="Gain_Loss"]<-"Gain+Loss"
plt_input$genomic_bin<-factor(plt_input$genomic_bin, levels=rev(c("Stable","Gain","Loss","Gain+Loss")))

names<-rev(c("Stable","Gain","Loss","Gain + Loss"))
col = c(Stable = "#EE3A2D", Gain = "#2270B6", Loss="#70ACD4", `Gain+Loss`="#193D6C")

plot<-ggplot(plt_input, aes(genomic_bin, fill=genomic_bin)) +
  geom_bar(aes(y = (..count..)/sum(..count..))) +
  #geom_text(aes(label=count), vjust=-0.2, size=3) +
  theme_bw() +
  scale_x_discrete(labels=c("Gain+Loss","Loss","Gain","Stable")) +
  coord_flip() +
  scale_fill_manual(values=col) +
  theme(axis.text.x = element_text(angle = 0, hjust=0, size=10),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        panel.border = element_rect(size=1),
        axis.text = element_text(color = "black")) +
  xlab(NULL) +
  ylab("Fraction Of Cases In Each Clonal Group") 

pdf("outputs/figure_2/clonal_group_fractions.pdf", width = 5, height = 3)
print(plot)
dev.off()


# read ATAC count matrix
count_matrix<-readRDS("inputs/bulk_count_matrix_raw.rds")
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

# Get PCA values to use in ggplot
pca_data<-plotPCA(rld, intgroup=c("timepoint", "patient"), returnData=TRUE)

ggplot(pca_data, aes(x=PC1,y=PC2,group=patient, color=patient, shape=timepoint)) +
  geom_point(size=3) +
  #scale_color_manual(values=c("darkblue","darkred","goldenrod")) +
  theme_classic()

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
deseq_Results<-as.data.frame(deseq_Results)
saveRDS(deseq_Results, "outputs/bulk_analysis/REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")
openxlsx::write.xlsx(rownames_to_column(deseq_Results),"outputs/bulk_analysis/REL_vs_DX_peak_accessibility_deseq_res_stable_cases.xlsx",row.names=F)

deseq_Results$color<-rep("stable",nrow(deseq_Results))
deseq_Results$color[deseq_Results$padj<0.05&deseq_Results$log2FoldChange>0.25]<-"up"
deseq_Results$color[deseq_Results$padj<0.05&deseq_Results$log2FoldChange<(-0.25)]<-"down"
deseq_Results$size<-as.numeric(sapply(strsplit(rownames(deseq_Results), "_"),`[`,3))-as.numeric(sapply(strsplit(rownames(deseq_Results), "_"),`[`,2))
deseq_Results<-deseq_Results[complete.cases(deseq_Results),]

pal<-c(stable="black",up="darkgreen",down="lightgreen")

plt<-ggplot(deseq_Results, aes(x=log2FoldChange,y=-log10(padj), color=color)) +
  ggrastr::rasterise(geom_point(size=0.5), dpi=500) +
  scale_color_manual(values=pal) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed", col="grey80") +
  geom_vline(xintercept = -0.25, linetype = "dashed", col="grey80") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col="grey80") +
  xlab("Log2FC REL - DX") +
  ylab("-log10(Adjusted p-val)") +
  ggtitle("Relapse vs Diagnosis Peaks Differential Accessibility") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color="black"))

pdf("outputs/figure_2/bulk_stable_rel_vs_dx_peaks_volcano.pdf", width = 4, height = 4)
print(plt)
dev.off()

# Heatmap
count_matrix_normalized<-normalize_count_matrix(count_matrix)
top_regions<-deseq_Results[deseq_Results$log2FoldChange>0&deseq_Results$padj<0.05,]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:200]
bottom_regions<-deseq_Results[deseq_Results$log2FoldChange<0&deseq_Results$padj<0.05,]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:200]
heatmap_input<-count_matrix_normalized[c(top_regions,bottom_regions),
                                       c(unique(bulk_sample_info_sub$name[bulk_sample_info_sub$timepoint=="DX"]),
                                         unique(bulk_sample_info_sub$name[bulk_sample_info_sub$timepoint=="REL"]))]

for(p in unique(bulk_sample_info_sub$patient)){
  message(p)
  heatmap_input[,unique(bulk_sample_info_sub$name[bulk_sample_info_sub$patient==p])]<-
    as.data.frame(t(scale(t(heatmap_input[,unique(bulk_sample_info_sub$name[bulk_sample_info_sub$patient==p])]))))
}

heatmap_input<-collapse_count_matrix(heatmap_input, bulk_sample_info_sub)
heatmap_input<-heatmap_input[c(top_regions,bottom_regions),
                             c(unique(bulk_sample_info_sub$sample[bulk_sample_info_sub$timepoint=="DX"]),
                               unique(bulk_sample_info_sub$sample[bulk_sample_info_sub$timepoint=="REL"]))]


if(T){
  annotations<-data.frame(bulk_sample_info_sub,stringsAsFactors = TRUE)
  annotations<-unique(annotations[,c("sample","timepoint")])
  rownames(annotations)<-NULL
  annotations<-column_to_rownames(annotations, "sample")
  annotations<-annotations[colnames(heatmap_input),,drop=FALSE]
  colnames(annotations)<-c("Timepoint")
  annotations$Timepoint<-as.factor(annotations$Timepoint)
  samples<-rownames(annotations)
  annotations<-as.data.frame(unclass(annotations))
  rownames(annotations)<-samples
  genotype_colors<-c("darkred","darkblue")
  names(genotype_colors)<-levels(annotations$Timepoint)
  anno_colors<-list(Timepoint=genotype_colors)
}

hm<-Heatmap(heatmap_input, 
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col=colorRamp2(seq(-0.75,0.75,length.out = 10), viridis(10)),
        show_row_names=FALSE, 
        show_column_names=T,
        heatmap_legend_param = list(title = "Row Z-Scores"),
        column_title="Stable Samples Differential Peak Accessibility",
        top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))

pdf("outputs/figure_2/bulk_stable_rel_vs_dx_peaks_heatmap.pdf", width = 5, height = 4)
print(hm)
dev.off()



## Gene Scores
gene_scores<-readRDS("inputs/bulk_gene_scores.rds")
gene_scores<-column_to_rownames(gene_scores,"Gene")
gene_scores<-floor(gene_scores)

# paired DEseq analysis - first do stable samples
bulk_sample_info_sub<-bulk_sample_info[bulk_sample_info$genotype_groups=="Stable",]
bulk_sample_info_sub<-bulk_sample_info_sub[!is.na(bulk_sample_info_sub$name),]
deseq_matrix<-gene_scores[,bulk_sample_info_sub$name]

data_reference<-bulk_sample_info_sub[,c("name","patient","timepoint")]
colnames(data_reference)<-c("sample","patient","timepoint")
dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~timepoint+patient))

# pca
rld <- varianceStabilizingTransformation(dds)
# Get PCA values to use in ggplot
pca_data<-plotPCA(rld, intgroup=c("timepoint", "patient"), returnData=TRUE)
pca_data$batch<-bulk_sample_info_sub$Batch[match(paste(pca_data$patient, pca_data$timepoint), paste(bulk_sample_info_sub$patient, bulk_sample_info_sub$timepoint))]

ggplot(pca_data, aes(x=PC1,y=PC2,group=patient, color=patient, shape=batch)) +
  geom_point(size=3) +
  #scale_color_manual(values=c("darkblue","darkred","goldenrod")) +
  theme_classic()

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
deseq_Results<-as.data.frame(deseq_Results)
saveRDS(deseq_Results, "outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")
openxlsx::write.xlsx(rownames_to_column(deseq_Results),"outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_stable_cases.xlsx",row.names=F)

# Deseq analysis of non-stable samples
bulk_sample_info_sub<-bulk_sample_info[bulk_sample_info$genotype_groups!="Stable",]
bulk_sample_info_sub<-bulk_sample_info_sub[!is.na(bulk_sample_info_sub$name),]
deseq_matrix<-gene_scores[,bulk_sample_info_sub$name]

data_reference<-bulk_sample_info_sub[,c("name","patient","timepoint")]
colnames(data_reference)<-c("sample","patient","timepoint")
dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~timepoint+patient))

# pca
rld <- varianceStabilizingTransformation(dds)
# Get PCA values to use in ggplot
pca_data<-plotPCA(rld, intgroup=c("timepoint", "patient"), returnData=TRUE)
pca_data$batch<-bulk_sample_info_sub$Batch[match(paste(pca_data$patient, pca_data$timepoint), paste(bulk_sample_info_sub$patient, bulk_sample_info_sub$timepoint))]

ggplot(pca_data, aes(x=PC1,y=PC2,group=patient, color=patient, shape=batch)) +
  geom_point(size=3) +
  #scale_color_manual(values=c("darkblue","darkred","goldenrod")) +
  theme_classic()

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
deseq_Results<-as.data.frame(deseq_Results)
saveRDS(deseq_Results, "outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_nonstable_cases.rds")
openxlsx::write.xlsx(rownames_to_column(deseq_Results),"outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_nonstable_cases.xlsx",row.names=F)

# plot differential gene accessibility for stable cases
deseq_Results<-readRDS("outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")

deseq_Results$color<-rep("stable",nrow(deseq_Results))
deseq_Results$color[deseq_Results$padj<0.05&deseq_Results$log2FoldChange>0.25]<-"up"
deseq_Results$color[deseq_Results$padj<0.05&deseq_Results$log2FoldChange<(-0.25)]<-"down"
deseq_Results<-deseq_Results[complete.cases(deseq_Results),]

pal<-c(stable="black",up="darkgreen",down="lightgreen")

plt<-ggplot(deseq_Results, aes(x=log2FoldChange,y=-log10(padj), color=color)) +
  ggrastr::rasterise(geom_point(size=0.5), dpi=500) +
  scale_color_manual(values=pal) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed", col="grey80") +
  geom_vline(xintercept = -0.25, linetype = "dashed", col="grey80") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col="grey80") +
  xlab("Log2FC REL - DX") +
  ylab("-log10(Adjusted p-val)") +
  ggtitle("Relapse vs Diagnosis Gene Differential Accessibility") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color="black"))

pdf("outputs/figure_2/bulk_stable_rel_vs_dx_gene_accessibility_volcano.pdf", width = 4, height = 4)
print(plt)
dev.off()


# Heatmap of differential stable genes in stable samples
bulk_sample_info_sub<-bulk_sample_info[bulk_sample_info$genotype_groups=="Stable",]
bulk_sample_info_sub<-bulk_sample_info_sub[!is.na(bulk_sample_info_sub$name),]

gene_scores<-readRDS("inputs/bulk_gene_scores.rds")
gene_scores<-column_to_rownames(gene_scores,"Gene")
count_matrix_normalized<-normalize_count_matrix(gene_scores)
top_regions<-deseq_Results[deseq_Results$log2FoldChange>0&deseq_Results$padj<0.05,]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:50]
bottom_regions<-deseq_Results[deseq_Results$log2FoldChange<0&deseq_Results$padj<0.05,]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:50]
heatmap_input<-count_matrix_normalized[c(top_regions,bottom_regions),
                                       c(unique(bulk_sample_info_sub$name[bulk_sample_info_sub$timepoint=="DX"]),
                                         unique(bulk_sample_info_sub$name[bulk_sample_info_sub$timepoint=="REL"]))]

for(p in unique(bulk_sample_info_sub$patient)){
  message(p)
  heatmap_input[,unique(bulk_sample_info_sub$name[bulk_sample_info_sub$patient==p])]<-
    as.data.frame(t(scale(t(heatmap_input[,unique(bulk_sample_info_sub$name[bulk_sample_info_sub$patient==p])]))))
}

heatmap_input<-collapse_count_matrix(heatmap_input, bulk_sample_info_sub)
heatmap_input<-heatmap_input[c(top_regions,bottom_regions),
                             c(unique(bulk_sample_info_sub$sample[bulk_sample_info_sub$timepoint=="DX"]),
                               unique(bulk_sample_info_sub$sample[bulk_sample_info_sub$timepoint=="REL"]))]


if(T){
  annotations<-data.frame(bulk_sample_info_sub,stringsAsFactors = TRUE)
  annotations<-unique(annotations[,c("sample","timepoint")])
  rownames(annotations)<-NULL
  annotations<-column_to_rownames(annotations, "sample")
  annotations<-annotations[colnames(heatmap_input),,drop=FALSE]
  colnames(annotations)<-c("Timepoint")
  annotations$Timepoint<-as.factor(annotations$Timepoint)
  samples<-rownames(annotations)
  annotations<-as.data.frame(unclass(annotations))
  rownames(annotations)<-samples
  genotype_colors<-c("darkred","darkblue")
  names(genotype_colors)<-levels(annotations$Timepoint)
  anno_colors<-list(Timepoint=genotype_colors)
}

hm<-Heatmap(heatmap_input, 
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col=colorRamp2(seq(-0.75,0.75,length.out = 10), viridis(10)),
        show_row_names=F, 
        show_column_names=T,
        heatmap_legend_param = list(title = "Row Z-Scores"),
        #column_title="Non-Stable Samples Differential Gene Accessibility",
        top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        left_annotation = rowAnnotation(genes = anno_mark(at = 1:20, 
                                                          labels = rownames(heatmap_input)[1:20], 
                                                          side="left",
                                                          labels_gp = gpar(fontsize = 8),
                                                          link_width =unit(5, "mm"),
                                                          padding = 0.7)),
        right_annotation = rowAnnotation(genes = anno_mark(at = 81:100, 
                                                           labels = rownames(heatmap_input)[81:100], 
                                                           side="right",
                                                           labels_gp = gpar(fontsize = 8),
                                                           link_width =unit(5, "mm"),
                                                           padding = 0.7)))

pdf("outputs/figure_2/bulk_stable_rel_vs_dx_gene_accessibility_heatmap.pdf", width = 6, height = 4)
print(hm)
dev.off()



# Heatmap of differential stable genes in nonstable samples
bulk_sample_info_sub<-bulk_sample_info[bulk_sample_info$genotype_groups!="Stable",]
bulk_sample_info_sub<-bulk_sample_info_sub[!is.na(bulk_sample_info_sub$name),]

gene_scores<-readRDS("inputs/bulk_gene_scores.rds")
gene_scores<-column_to_rownames(gene_scores,"Gene")
count_matrix_normalized<-normalize_count_matrix(gene_scores)
top_regions<-deseq_Results[deseq_Results$log2FoldChange>0&deseq_Results$padj<0.05,]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:50]
bottom_regions<-deseq_Results[deseq_Results$log2FoldChange<0&deseq_Results$padj<0.05,]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:50]
heatmap_input<-count_matrix_normalized[c(top_regions,bottom_regions),
                                       c(unique(bulk_sample_info_sub$name[bulk_sample_info_sub$timepoint=="DX"]),
                                         unique(bulk_sample_info_sub$name[bulk_sample_info_sub$timepoint=="REL"]))]

for(p in unique(bulk_sample_info_sub$patient)){
  message(p)
  heatmap_input[,unique(bulk_sample_info_sub$name[bulk_sample_info_sub$patient==p])]<-
    as.data.frame(t(scale(t(heatmap_input[,unique(bulk_sample_info_sub$name[bulk_sample_info_sub$patient==p])]))))
}

heatmap_input<-collapse_count_matrix(heatmap_input, bulk_sample_info_sub)
heatmap_input<-heatmap_input[c(top_regions,bottom_regions),
                             c(unique(bulk_sample_info_sub$sample[bulk_sample_info_sub$timepoint=="DX"]),
                               unique(bulk_sample_info_sub$sample[bulk_sample_info_sub$timepoint=="REL"]))]


if(T){
  annotations<-data.frame(bulk_sample_info_sub,stringsAsFactors = TRUE)
  annotations<-unique(annotations[,c("sample","timepoint")])
  rownames(annotations)<-NULL
  annotations<-column_to_rownames(annotations, "sample")
  annotations<-annotations[colnames(heatmap_input),,drop=FALSE]
  colnames(annotations)<-c("Timepoint")
  annotations$Timepoint<-as.factor(annotations$Timepoint)
  samples<-rownames(annotations)
  annotations<-as.data.frame(unclass(annotations))
  rownames(annotations)<-samples
  genotype_colors<-c("darkred","darkblue")
  names(genotype_colors)<-levels(annotations$Timepoint)
  anno_colors<-list(Timepoint=genotype_colors)
}

hm<-Heatmap(heatmap_input, 
            cluster_columns = FALSE,
            cluster_rows = FALSE,
            col=colorRamp2(seq(-0.75,0.75,length.out = 10), viridis(10)),
            show_row_names=F, 
            show_column_names=T,
            heatmap_legend_param = list(title = "Row Z-Scores"),
            #column_title="Non-Stable Samples Differential Gene Accessibility",
            top_annotation=HeatmapAnnotation(df=annotations, col=anno_colors, annotation_name_side = "left"),
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8),
            left_annotation = rowAnnotation(genes = anno_mark(at = 1:20, 
                                                              labels = rownames(heatmap_input)[1:20], 
                                                              side="left",
                                                              labels_gp = gpar(fontsize = 8),
                                                              link_width =unit(5, "mm"),
                                                              padding = 0.7)),
            right_annotation = rowAnnotation(genes = anno_mark(at = 81:100, 
                                                               labels = rownames(heatmap_input)[81:100], 
                                                               side="right",
                                                               labels_gp = gpar(fontsize = 8),
                                                               link_width =unit(5, "mm"),
                                                               padding = 0.7)))

pdf("outputs/figure_2/bulk_non_stable_rel_vs_dx_gene_accessibility_heatmap_diff_stable_genes.pdf", width = 6, height = 4)
print(hm)
dev.off()




## differential analysis of bulk separate patients for downstream single cell analysis
## Gene Scores
gene_scores<-readRDS("inputs/bulk_gene_scores.rds")
gene_scores<-column_to_rownames(gene_scores,"Gene")
gene_scores<-floor(gene_scores)

for(p in c("SU142", "SU360", "SU892")){
  
  # paired DEseq analysis - first do stable samples
  bulk_sample_info_sub<-bulk_sample_info[bulk_sample_info$genotype_groups=="Stable"&bulk_sample_info$patient==p,]
  bulk_sample_info_sub<-bulk_sample_info_sub[!is.na(bulk_sample_info_sub$name),]
  deseq_matrix<-gene_scores[,bulk_sample_info_sub$name]
  
  data_reference<-bulk_sample_info_sub[,c("name","timepoint")]
  colnames(data_reference)<-c("sample","timepoint")
  dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~timepoint))
  
  # pca
  rld <- varianceStabilizingTransformation(dds)
  # Get PCA values to use in ggplot
  # pca_data<-plotPCA(rld, intgroup=c("timepoint"), returnData=TRUE)
  # pca_data$batch<-bulk_sample_info_sub$Batch[match(paste(pca_data$patient, pca_data$timepoint), paste(bulk_sample_info_sub$patient, bulk_sample_info_sub$timepoint))]
  # 
  # ggplot(pca_data, aes(x=PC1,y=PC2,group=patient, color=patient, shape=batch)) +
  #   geom_point(size=3) +
  #   #scale_color_manual(values=c("darkblue","darkred","goldenrod")) +
  #   theme_classic()
  
  dds<-estimateSizeFactors(dds)
  dds <- DESeq(dds)
  deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
  deseq_Results<-as.data.frame(deseq_Results)
  saveRDS(deseq_Results, paste0("outputs/bulk_analysis/",p,"_REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds"))
  write.xlsx(rownames_to_column(deseq_Results),paste0("outputs/bulk_analysis/",p,"_REL_vs_DX_gene_accessibility_deseq_res_stable_cases.xlsx"),row.names=F)
}
