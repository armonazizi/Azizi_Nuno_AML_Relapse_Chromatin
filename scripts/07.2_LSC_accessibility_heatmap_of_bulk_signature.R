# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 3
# LSC peak fold change vs Bulk/Blast peak fold change

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")

lsc_cases<-unique(sample_info$patient[sample_info$type=="LSC"])
lsc_sample_info<-sample_info[sample_info$patient%in%lsc_cases&sample_info$type%in%c("LSC","Blast","Bulk"),]

genotype_groups<-genotyping_info[,c("ID","genomic_bin")]
lsc_sample_info$genotype_groups<-genotype_groups$genomic_bin[match(lsc_sample_info$patient,genotype_groups$ID)]

lsc_sample_info_sub<-lsc_sample_info[lsc_sample_info$genotype_groups=="Stable"&lsc_sample_info$type=="LSC",]
lsc_sample_info_sub<-lsc_sample_info_sub[!is.na(lsc_sample_info_sub$name),]

bulk_sample_info<-sample_info[sample_info$type%in%c("LSC","Blast","Bulk"),]

genotype_groups<-genotyping_info[,c("ID","genomic_bin")]
bulk_sample_info$genotype_groups<-genotype_groups$genomic_bin[match(bulk_sample_info$patient,genotype_groups$ID)]

bulk_sample_info_sub<-bulk_sample_info[bulk_sample_info$genotype_groups=="Stable"&bulk_sample_info$type=="Bulk",]
bulk_sample_info_sub<-bulk_sample_info_sub[!is.na(bulk_sample_info_sub$name),]


# Heatmap of differential stable genes in stable samples
deseq_Results<-readRDS("outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")
deseq_Results<-deseq_Results[complete.cases(deseq_Results),]

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

#pdf("outputs/figure_2/bulk_stable_rel_vs_dx_gene_accessibility_heatmap.pdf", width = 6, height = 4)
print(hm)
#dev.off()



# Heatmap of differential stable genes in LSC stable samples
gene_scores<-readRDS("inputs/bulk_gene_scores.rds")
gene_scores<-column_to_rownames(gene_scores,"Gene")
count_matrix_normalized<-normalize_count_matrix(gene_scores)
top_regions<-deseq_Results[deseq_Results$log2FoldChange>0&deseq_Results$padj<0.05,]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:50]
bottom_regions<-deseq_Results[deseq_Results$log2FoldChange<0&deseq_Results$padj<0.05,]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:50]
heatmap_input<-count_matrix_normalized[c(top_regions,bottom_regions),
                                       c(unique(lsc_sample_info_sub$name[lsc_sample_info_sub$timepoint=="DX"]),
                                         unique(lsc_sample_info_sub$name[lsc_sample_info_sub$timepoint=="REL"]))]

for(p in unique(lsc_sample_info_sub$patient)){
  message(p)
  heatmap_input[,unique(lsc_sample_info_sub$name[lsc_sample_info_sub$patient==p])]<-
    as.data.frame(t(scale(t(heatmap_input[,unique(lsc_sample_info_sub$name[lsc_sample_info_sub$patient==p])]))))
}

heatmap_input<-collapse_count_matrix(heatmap_input, lsc_sample_info_sub)
heatmap_input<-heatmap_input[c(top_regions,bottom_regions),
                             c(unique(lsc_sample_info_sub$sample[lsc_sample_info_sub$timepoint=="DX"]),
                               unique(lsc_sample_info_sub$sample[lsc_sample_info_sub$timepoint=="REL"]))]


if(T){
  annotations<-data.frame(lsc_sample_info_sub,stringsAsFactors = TRUE)
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

pdf("outputs/figure_2/lsc_stable_rel_vs_dx_gene_accessibility_heatmap_diff_stable_genes.pdf", width = 6, height = 4)
print(hm)
dev.off()