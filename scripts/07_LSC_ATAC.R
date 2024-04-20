# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 3
# Chromatin accessibility analysis of LSCs and Blasts sorted from diagnosis and relapse AML


setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

# Read metadata and genotyping info
sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/genotyping_formatted_V2.xlsx")

lsc_cases<-unique(sample_info$patient[sample_info$type=="LSC"])
lsc_sample_info<-sample_info[sample_info$patient%in%lsc_cases&sample_info$type%in%c("LSC","Blast"),]

genotype_groups<-genotyping_info[,c("ID","genomic_bin")]
lsc_sample_info$genotype_groups<-genotype_groups$genomic_bin[match(lsc_sample_info$patient,genotype_groups$ID)]

# Compare number of cases with LSCs in stable and nonstable groups
sample_info$genotype_groups<-genotype_groups$genomic_bin[match(sample_info$patient,genotype_groups$ID)]
lsc_case_counts<-sample_info[,c("patient","type", "genotype_groups")]
lsc_pt<-unique(lsc_case_counts$patient[lsc_case_counts$type=="LSC"])
non_lsc_pt<-setdiff(unique(lsc_case_counts$patient),lsc_pt)
lsc_case_counts<-data.frame(patient=c(lsc_pt,non_lsc_pt),
                    has_lsc=c(rep("LSCs Present",length(lsc_pt)),rep("LSCs Absent",length(non_lsc_pt))))
lsc_case_counts$genotype_groups<-genotype_groups$genomic_bin[match(lsc_case_counts$patient,genotype_groups$ID)]
lsc_case_counts$genotype_groups[lsc_case_counts$genotype_groups!="Stable"]<-"Non-Stable"
plt_input<-reshape2::melt(as.matrix(table(lsc_case_counts[,c("has_lsc", "genotype_groups")])))
plt<-ggplot(plt_input, aes(x=genotype_groups, y=value, fill=has_lsc)) +
  geom_bar(stat="identity") +
  theme_classic() +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_fill_manual(values = c("grey30","grey70")) +
  ggtitle("Number Of Cases With A Detectable LSC Population") +
  ylab("Number Of Cases") +
  annotate("text", x=2, y=12, label= "p = 0.0152") +
  theme(legend.title = element_blank())

pdf("outputs/figure_3_LSC/stable_vs_nonstable_lsc_case_counts.pdf", width = 5, height = 4)
plt
dev.off()



# read ATAC count matrix
count_matrix<-readRDS("inputs/bulk_count_matrix_raw.rds")
rownames(count_matrix)<-paste(count_matrix$Chr,count_matrix$Start,count_matrix$End,sep="_")
count_matrix<-count_matrix[,4:length(colnames(count_matrix))]

count_matrix_normalized<-normalize_count_matrix(count_matrix)
count_matrix_normalized<-count_matrix_normalized[,lsc_sample_info$name]

# calculate dx-rel similarity
lsc_sample_info$sample_type<-paste(lsc_sample_info$sample, lsc_sample_info$type, sep="_")
lsc_sample_info$case_type<-paste(lsc_sample_info$patient, lsc_sample_info$type, sep="_")
lsc_sample_info$case_timepoint<-paste(lsc_sample_info$patient, lsc_sample_info$timepoint, sep="_")

count_matrix_collapsed<-collapse_count_matrix(count_matrix_normalized, lsc_sample_info, "sample_type")

correlations<-cor(count_matrix_collapsed)
corr_df<-data.frame()

for(p in unique(lsc_sample_info$case_type)){
  message(p)
  corr<-correlations[unique(lsc_sample_info$sample_type[lsc_sample_info$case_type==p&lsc_sample_info$timepoint=="DX"][1]),
                     unique(lsc_sample_info$sample_type[lsc_sample_info$case_type==p&lsc_sample_info$timepoint=="REL"][1])]
  corr_df[p,"Correlation"]<-corr
}
corr_df$genotype_group<-lsc_sample_info$genotype_groups[match(rownames(corr_df), lsc_sample_info$case_type)]
corr_df<-corr_df[complete.cases(corr_df),]
corr_df$type<-sapply(strsplit(rownames(corr_df), "_"),'[',2)
corr_df$patient<-sapply(strsplit(rownames(corr_df), "_"),'[',1)

plt<-ggplot(corr_df,aes(x=type, y=Correlation,group=patient, color=genotype_group)) +
  geom_boxplot(mapping = aes(x=type, y=Correlation), inherit.aes = F, size=1) +
  geom_point() +
  geom_line(size=1) +
  scale_color_manual(values=c("darkred","darkblue","goldenrod")) +
  xlab("Cell Type") +
  ylab("DX vs REL ATAC Similarity") +
  theme_classic() +
  stat_compare_means(method="t.test", paired = TRUE, comparisons = list(c("LSC","Blast")))

print(plt)

pdf("outputs/figure_3_LSC/LSC_vs_Blast_relapse_accessibility_change.pdf", width = 5, height = 4)
plt
dev.off()



#### compare LSCs vs blast similarity at diagnosis vs relapse ####
correlations<-cor(count_matrix_collapsed)
corr_df<-data.frame()

for(p in unique(lsc_sample_info$case_timepoint)){
  message(p)
  corr<-correlations[unique(lsc_sample_info$sample_type[lsc_sample_info$case_timepoint==p&lsc_sample_info$type=="Blast"][1]),
                     unique(lsc_sample_info$sample_type[lsc_sample_info$case_timepoint==p&lsc_sample_info$type=="LSC"][1])]
  corr_df[p,"Correlation"]<-corr
}
corr_df$genotype_group<-lsc_sample_info$genotype_groups[match(rownames(corr_df), lsc_sample_info$case_timepoint)]
corr_df<-corr_df[complete.cases(corr_df),]
corr_df$type<-sapply(strsplit(rownames(corr_df), "_"),'[',2)
corr_df$patient<-sapply(strsplit(rownames(corr_df), "_"),'[',1)

corr_df$genotype_group[corr_df$genotype_group!="Stable"]<-"Unstable"

plt<-ggplot(corr_df,aes(x=type, y=Correlation,group=patient, color=genotype_group)) +
  #geom_boxplot(mapping = aes(x=type, y=Correlation), inherit.aes = F, size=1) +
  geom_line(size=0.5) +
  geom_point(shape=21, fill="white") +
  scale_color_manual(name = "Genetic\nEvolution", values=c("red","darkblue")) +
  xlab("Timepoint") +
  ylab("LSC vs Non-LSC ATAC Similarity") +
  theme_classic() +
  stat_compare_means(method="t.test", paired = TRUE, comparisons = list(c("DX","REL")))

print(plt)

pdf("outputs/figure_3_LSC/DX_vs_REL_LSC_vs_Blast_similarity.pdf", width = 3, height = 3.5)
plt
dev.off()






# Do LSC REL vs DX differential accessibility

# peaks
# read ATAC count matrix
count_matrix<-readRDS("inputs/bulk_count_matrix_raw.rds")
rownames(count_matrix)<-paste(count_matrix$Chr,count_matrix$Start,count_matrix$End,sep="_")
count_matrix<-count_matrix[,4:length(colnames(count_matrix))]

## paired DEseq analysis of stable samples
lsc_sample_info_sub<-lsc_sample_info[lsc_sample_info$genotype_groups=="Stable",]
lsc_sample_info_sub<-lsc_sample_info_sub[!is.na(lsc_sample_info_sub$name),]
deseq_matrix<-count_matrix[,lsc_sample_info_sub$name]

data_reference<-lsc_sample_info_sub[,c("name","patient","timepoint")]
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
saveRDS(deseq_Results, "outputs/lsc_analysis/LSC_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")



## Gene Scores
gene_scores<-readRDS("inputs/bulk_gene_scores.rds")
gene_scores<-column_to_rownames(gene_scores,"Gene")
gene_scores<-floor(gene_scores)

# paired DEseq analysis - first do stable samples
deseq_matrix<-gene_scores[,lsc_sample_info_sub$name]

data_reference<-lsc_sample_info_sub[,c("name","patient","timepoint")]
colnames(data_reference)<-c("sample","patient","timepoint")
dds<-DESeqDataSetFromMatrix(deseq_matrix,data_reference, formula(~timepoint+patient))

# pca
rld <- varianceStabilizingTransformation(dds)
# Get PCA values to use in ggplot
pca_data<-plotPCA(rld, intgroup=c("timepoint", "patient"), returnData=TRUE)
pca_data$batch<-lsc_sample_info_sub$Batch[match(paste(pca_data$patient, pca_data$timepoint), paste(lsc_sample_info_sub$patient, lsc_sample_info_sub$timepoint))]

ggplot(pca_data, aes(x=PC1,y=PC2,group=patient, color=patient, shape=batch)) +
  geom_point(size=3) +
  #scale_color_manual(values=c("darkblue","darkred","goldenrod")) +
  theme_classic()

dds<-estimateSizeFactors(dds)
dds <- DESeq(dds)
deseq_Results <- results(dds, contrast=c("timepoint","REL","DX"))
deseq_Results<-as.data.frame(deseq_Results)
saveRDS(deseq_Results, "outputs/lsc_analysis/LSC_REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")



## Heatmap of differential gene accessibility for LSCs stable cases
deseq_Results<-readRDS("outputs/lsc_analysis/LSC_REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")

lsc_sample_info_sub<-lsc_sample_info[lsc_sample_info$genotype_groups=="Stable",]
lsc_sample_info_sub<-lsc_sample_info_sub[!is.na(lsc_sample_info_sub$name),]

gene_scores<-readRDS("inputs/bulk_gene_scores.rds")
gene_scores<-column_to_rownames(gene_scores,"Gene")
count_matrix_normalized<-normalize_count_matrix(gene_scores)
top_regions<-deseq_Results[deseq_Results$log2FoldChange>0&deseq_Results$padj<0.05,]
top_regions<-top_regions[complete.cases(top_regions),]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))]#[1:50]
bottom_regions<-deseq_Results[deseq_Results$log2FoldChange<0&deseq_Results$padj<0.05,]
bottom_regions<-bottom_regions[complete.cases(bottom_regions),]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)]#[1:50]
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

print(hm)

# pdf("outputs/figure_2/bulk_stable_rel_vs_dx_gene_accessibility_heatmap.pdf", width = 6, height = 4)
# print(hm)
# dev.off()






# Do volcano plot of LSC differential analysis and map differential bulk gene signature

# peaks
bulk_diff_res<-readRDS("outputs/bulk_analysis/REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")
lsc_diff_res<-readRDS("outputs/lsc_analysis/LSC_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")

bulk_diff_res<-bulk_diff_res[complete.cases(bulk_diff_res),]
lsc_diff_res<-lsc_diff_res[complete.cases(lsc_diff_res),]

top_regions<-bulk_diff_res[bulk_diff_res$log2FoldChange>0&bulk_diff_res$padj<0.05,]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:200]
bottom_regions<-bulk_diff_res[bulk_diff_res$log2FoldChange<0&bulk_diff_res$padj<0.05,]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:200]

bulk_diff_regions<-c(top_regions, bottom_regions)

top_regions<-lsc_diff_res[lsc_diff_res$log2FoldChange>0&lsc_diff_res$padj<0.05,]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:200]
bottom_regions<-lsc_diff_res[lsc_diff_res$log2FoldChange<0&lsc_diff_res$padj<0.05,]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:200]

lsc_diff_regions<-c(top_regions, bottom_regions)


lsc_diff_res<-lsc_diff_res[complete.cases(lsc_diff_res),]

bulk_diff_res$color<-rep("",nrow(bulk_diff_res))
bulk_diff_res$color[rownames(bulk_diff_res)%in%lsc_diff_regions]<-"LSC Differential"
bulk_diff_res$color[rownames(bulk_diff_res)%in%bulk_diff_regions]<-"Bulk Differential"

plt<-ggplot(bulk_diff_res, aes(x=log2FoldChange,y=-log10(padj), color=color)) +
  ggrastr::rasterise(geom_point(size=0.5, alpha=0.5), dpi=300) +
  geom_point(data = bulk_diff_res[lsc_diff_regions,], size=0.5, color="red") +
  geom_point(data = bulk_diff_res[bulk_diff_regions,], size=0.5, color="blue") +
  scale_color_manual(values=c("grey30","blue","red")) +
  theme_bw() +
  xlab("Log2FC REL - DX") +
  ylab("-log10(Adjusted p-val)") +
  ggtitle("Bulk Relapse vs Diagnosis Peak Acessibility") +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_line(colour = "black", size = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(),
        axis.text = element_text(color="black"))

print(plt)

pdf("outputs/figure_3_LSC/bulk_volcano_genes_mapped.pdf", width = 5, height = 4)
plt
dev.off()


lsc_diff_res$color<-rep("",nrow(lsc_diff_res))
lsc_diff_res$color[rownames(lsc_diff_res)%in%lsc_diff_regions]<-"LSC Differential"
lsc_diff_res$color[rownames(lsc_diff_res)%in%bulk_diff_regions]<-"Bulk Differential"

plt<-ggplot(lsc_diff_res, aes(x=log2FoldChange,y=-log10(padj), color=color)) +
  ggrastr::rasterise(geom_point(size=0.5, alpha=0.5), dpi=300) +
  geom_point(data = lsc_diff_res[lsc_diff_regions,], size=0.5, color="red") +
  geom_point(data = lsc_diff_res[bulk_diff_regions,], size=0.5, color="blue") +
  scale_color_manual(values=c("grey30","black","darkred")) +
  theme_bw() +
  xlab("Log2FC REL - DX") +
  ylab("-log10(Adjusted p-val)") +
  ggtitle("LSC Relapse vs Diagnosis Peak Acessibility") +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_line(colour = "black", size = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(),
        axis.text = element_text(color="black"))

print(plt)

pdf("outputs/figure_3_LSC/lsc_volcano_genes_mapped.pdf", width = 5, height = 4)
plt
dev.off()





# Do volcano plot of LSC differential analysis and map differential bulk gene signature

# peaks
blast_diff_res<-readRDS("outputs/lsc_analysis/Blast_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")
lsc_diff_res<-readRDS("outputs/lsc_analysis/LSC_REL_vs_DX_peak_accessibility_deseq_res_stable_cases.rds")

blast_diff_res<-blast_diff_res[complete.cases(blast_diff_res),]
lsc_diff_res<-lsc_diff_res[complete.cases(lsc_diff_res),]

top_regions<-blast_diff_res[blast_diff_res$log2FoldChange>0&blast_diff_res$padj<0.05,]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:200]
bottom_regions<-blast_diff_res[blast_diff_res$log2FoldChange<0&blast_diff_res$padj<0.05,]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:200]

blast_diff_regions<-c(top_regions, bottom_regions)

top_regions<-lsc_diff_res[lsc_diff_res$log2FoldChange>0&lsc_diff_res$padj<0.05,]
top_regions<-rownames(top_regions)[rev(order(top_regions$log2FoldChange))][1:200]
bottom_regions<-lsc_diff_res[lsc_diff_res$log2FoldChange<0&lsc_diff_res$padj<0.05,]
bottom_regions<-rownames(bottom_regions)[order(bottom_regions$log2FoldChange)][1:200]

lsc_diff_regions<-c(top_regions, bottom_regions)


lsc_diff_res<-lsc_diff_res[complete.cases(lsc_diff_res),]

blast_diff_res$color<-rep("",nrow(blast_diff_res))
blast_diff_res$color[rownames(blast_diff_res)%in%lsc_diff_regions]<-"LSC Differential"
blast_diff_res$color[rownames(blast_diff_res)%in%blast_diff_regions]<-"Blast Differential"

plt<-ggplot(blast_diff_res, aes(x=log2FoldChange,y=-log10(padj), color=color)) +
  ggrastr::rasterise(geom_point(size=0.5, alpha=0.5), dpi=300) +
  geom_point(data = blast_diff_res[lsc_diff_regions,], size=0.5, color="red") +
  geom_point(data = blast_diff_res[blast_diff_regions,], size=0.5, color="blue") +
  scale_color_manual(values=c("grey30","blue","red")) +
  theme_bw() +
  xlab("Log2FC REL - DX") +
  ylab("-log10(Adjusted p-val)") +
  ggtitle("Blast Relapse vs Diagnosis Peak Acessibility") +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_line(colour = "black", size = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(),
        axis.text = element_text(color="black"))

print(plt)

pdf("outputs/figure_3_LSC/blast_volcano_lsc_genes_mapped.pdf", width = 5, height = 4)
plt
dev.off()


lsc_diff_res$color<-rep("",nrow(lsc_diff_res))
lsc_diff_res$color[rownames(lsc_diff_res)%in%lsc_diff_regions]<-"LSC Differential"
lsc_diff_res$color[rownames(lsc_diff_res)%in%blast_diff_regions]<-"Blast Differential"

plt<-ggplot(lsc_diff_res, aes(x=log2FoldChange,y=-log10(padj), color=color)) +
  ggrastr::rasterise(geom_point(size=0.5, alpha=0.5), dpi=300) +
  geom_point(data = lsc_diff_res[lsc_diff_regions,], size=0.5, color="red") +
  geom_point(data = lsc_diff_res[blast_diff_regions,], size=0.5, color="blue") +
  scale_color_manual(values=c("grey30","black","darkred")) +
  theme_bw() +
  xlab("Log2FC REL - DX") +
  ylab("-log10(Adjusted p-val)") +
  ggtitle("LSC Relapse vs Diagnosis Peak Acessibility") +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_line(colour = "black", size = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(),
        axis.text = element_text(color="black"))

print(plt)

pdf("outputs/figure_3_LSC/lsc_volcano_blast_genes_mapped.pdf", width = 5, height = 4)
plt
dev.off()












# Find fraction of LSC significant peaks that are also significant in blasts and vice versa
lsc_significant_bulk_significant<-length(intersect(rownames(lsc_diff_res[lsc_diff_res$padj<0.05,]),
                                                    rownames(bulk_diff_res[bulk_diff_res$padj<0.05,])))
lsc_significant_bulk_nonsignificant<-sum(lsc_diff_res$padj<0.05)-lsc_significant_bulk_significant

bulk_significant_lsc_significant<-length(intersect(rownames(bulk_diff_res[bulk_diff_res$padj<0.05,]),
                                                    rownames(lsc_diff_res[lsc_diff_res$padj<0.05,])))
bulk_significant_lsc_nonsignificant<-sum(bulk_diff_res$padj<0.05)-bulk_significant_lsc_significant

plt_input<-data.frame(CellType=c("LSC","LSC","Bulk","Bulk"),
                      Significance=c("Significant In Other Cell Type","Not Significant In Other Cell Type","Significant In Other Cell Type","Not Significant In Other Cell Type"),
                      significant_peaks=c(lsc_significant_bulk_significant,lsc_significant_bulk_nonsignificant,bulk_significant_lsc_significant,bulk_significant_lsc_nonsignificant))

plt<-ggplot(plt_input, aes(fill=Significance, y=significant_peaks, x=CellType)) + 
  geom_bar(position="fill", stat="identity", color="black", size=1) +
  scale_fill_manual(values = c("grey30","grey80")) +
  theme_classic() +
  ylab("Fraction Of Significant Peaks") +
  annotate("text", x=1.5, y=1.1, label= "****", size=5) +
  annotate("text", x=2, y=1.05, label= "n=460", size=3) +
  annotate("text", x=1, y=1.05, label= "n=11300", size=3)

print(plt)

pdf("outputs/figure_3_LSC/lsc_vs_bulk_significant_gene_overlaps.pdf", width = 5, height = 4)
plt
dev.off()

