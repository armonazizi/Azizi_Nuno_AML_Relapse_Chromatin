# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 2
# Summarization of cibersort fractions analysis of bulk AML relapse ATAC data

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

sample_info<-read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")
genotyping_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/genotyping_final.xlsx")
genotype_groups<-genotyping_info[,c("Sample.ID","Genomic.Bin")]

# read CSX fractions
cibersort_fractions<-as.data.frame(fread("inputs/cibersort_outputs/AML_vs_normal_heme_signature_CSX_fractions.txt"))

cibersort_fractions$`B/Plasma`<-rowSums(cibersort_fractions[,c("Bcell"),drop=F])
cibersort_fractions$CLP<-rowSums(cibersort_fractions[,c("CLP"),drop=F])
cibersort_fractions$`CMP/LMPP`<-rowSums(cibersort_fractions[,c("LMPP", "CMP")])
cibersort_fractions$Erythroid<-rowSums(cibersort_fractions[,c("Ery","MEP")])
cibersort_fractions$GMP<-rowSums(cibersort_fractions[,c("GMP"),drop=F])
cibersort_fractions$HSC<-rowSums(cibersort_fractions[,c("HSC"),drop=F])
cibersort_fractions$`Mono/DC/Baso`<-rowSums(cibersort_fractions[,c("Mono"),drop=F])
cibersort_fractions$`T/NK`<-rowSums(cibersort_fractions[,c("CD4Tcell","CD8Tcell","NKcell")])
cibersort_fractions$sample<-sample_info$sample[match(sample_info$name, cibersort_fractions$Mixture)]
cibersort_fractions$is_bulk<-sample_info$is_bulk[match(sample_info$name, cibersort_fractions$Mixture)]

classification_celltypes<-c("HSC","CMP/LMPP","CLP","GMP","B/Plasma","Erythroid","Mono/DC/Baso", "T/NK")


# fraction barplots for bulk samples
bulk_samples<-sample_info$ID[sample_info$is_bulk==1]
bulk_sample_info<-sample_info[sample_info$ID%in%bulk_samples,]
bulk_sample_info$genotype_groups<-genotype_groups$Genomic.Bin[match(bulk_sample_info$patient,genotype_groups$Sample.ID)]

bulk_fractions_sub<-column_to_rownames(cibersort_fractions, "Mixture")[bulk_sample_info$name, classification_celltypes]

bulk_fractions_sub$sample<-substr(rownames(bulk_fractions_sub),1,nchar(rownames(bulk_fractions_sub))-2)
bulk_fractions_sub<-aggregate(bulk_fractions_sub[,classification_celltypes], by=list(bulk_fractions_sub$sample), FUN=mean)
bulk_fractions_sub<-column_to_rownames(bulk_fractions_sub, "Group.1")

bulk_fractions_sub$patient<-strsplit(rownames(bulk_fractions_sub),"_") %>% sapply('[',1)
bulk_fractions_sub$timepoint<-strsplit(rownames(bulk_fractions_sub),"_") %>% sapply('[',2)
bulk_fractions_sub$type<-strsplit(rownames(bulk_fractions_sub),"_") %>% sapply('[',3)
bulk_fractions_sub$patient_timepoint_type<-paste(bulk_fractions_sub$patient,bulk_fractions_sub$timepoint,bulk_fractions_sub$type, sep="_")

bulk_fractions_sub_melt<-reshape2::melt(as.matrix(bulk_fractions_sub[,classification_celltypes]))
bulk_fractions_sub_melt$Var2<-paste0(bulk_fractions_sub_melt$Var2, "-like")

celltype_pal<-c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB")
names(celltype_pal)<-c("B/Plasma-like","CLP-like","CMP/LMPP-like","Erythroid-like","GMP-like","HSC-like","Mono/DC/Baso-like","T/NK-like")

plt<-ggplot(data=bulk_fractions_sub_melt, aes(x=Var1, y=value, fill=Var2)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(name="CSX Fractions",values = celltype_pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("") +
  ylab("Fractional Contribution")

pdf(paste0("outputs/figure_2/csx_fractions_individual_samples.pdf"), width = 10, height = 4)
print(plt)
dev.off()


# average all samples and replot
bulk_fractions_sub_avg<-aggregate(bulk_fractions_sub[,classification_celltypes], by=list(bulk_fractions_sub$timepoint), FUN=mean)
bulk_fractions_sub_avg<-column_to_rownames(bulk_fractions_sub_avg, "Group.1")

bulk_fractions_sub_avg_melt<-reshape2::melt(as.matrix(bulk_fractions_sub_avg[,classification_celltypes]))
bulk_fractions_sub_avg_melt$Var2<-paste0(bulk_fractions_sub_avg_melt$Var2, "-like")

celltype_pal<-c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB")
names(celltype_pal)<-c("B/Plasma-like","CLP-like","CMP/LMPP-like","Erythroid-like","GMP-like","HSC-like","Mono/DC/Baso-like","T/NK-like")

plt<-ggplot(data=bulk_fractions_sub_avg_melt, aes(x=Var1, y=value, fill=Var2)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(name="CSX Fractions",values = celltype_pal) +
  theme_classic() +
  xlab("") +
  ylab("Fractional Contribution")

pdf(paste0("outputs/figure_2/csx_fractions_all_samples_avg.pdf"), width = 4, height = 4)
print(plt)
dev.off()


# cibersort closest normal summary - all samples
cibersort_fractions$sample<-substr(cibersort_fractions$Mixture,1,nchar(cibersort_fractions$Mixture)-2)
cibersort_fractions<-aggregate(cibersort_fractions[,classification_celltypes], by=list(cibersort_fractions$sample), FUN=mean)

cibersort_closest_normal<-data.frame(mixture=cibersort_fractions$Group.1,
                                     closest_norm=classification_celltypes[max.col(cibersort_fractions[,classification_celltypes],ties.method="first")])


cibersort_closest_normal$patient<-strsplit(cibersort_closest_normal$mixture,"_") %>% sapply('[',1)
cibersort_closest_normal$timepoint<-strsplit(cibersort_closest_normal$mixture,"_") %>% sapply('[',2)
cibersort_closest_normal$type<-strsplit(cibersort_closest_normal$mixture,"_") %>% sapply('[',3)
cibersort_closest_normal$patient_type<-paste(cibersort_closest_normal$patient,cibersort_closest_normal$type, sep="_")
cibersort_closest_normal$closest_norm<-paste0(cibersort_closest_normal$closest_norm, "-like")

heatmap_input<-dcast(cibersort_closest_normal, patient_type ~ timepoint, value.var = "closest_norm")
heatmap_input<-column_to_rownames(heatmap_input, "patient_type")
heatmap_input<-heatmap_input[order(heatmap_input$REL),]
heatmap_input<-heatmap_input[order(heatmap_input$DX),]


annotations<-data.frame(type=strsplit(rownames(heatmap_input),"_") %>% sapply('[',2))
annotations$type<-tolower(annotations$type)

heatmap_input<-heatmap_input[order(annotations$type),]
annotations<-annotations[order(annotations$type),,drop=F]

celltype_pal<-c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB")
names(celltype_pal)<-c("B/Plasma-like","CLP-like","CMP/LMPP-like","Erythroid-like","GMP-like","HSC-like","Mono/DC/Baso-like","T/NK-like")

celltype_pal<-celltype_pal[names(celltype_pal)%in%unique(unlist(heatmap_input))]
anno_colors<-list(type=c(lsc="yellow",bulk="red",blast="darkblue"))

plt<-ComplexHeatmap::Heatmap(heatmap_input, 
                        col=celltype_pal,
                        show_row_names=T, 
                        show_column_names=T, 
                        cluster_rows = T,
                        cluster_columns = F,
                        heatmap_legend_param = list(title = "Number Overlaps"),
                        column_title="Cibersort Closest Normals",
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        rect_gp = gpar(col = "white", lwd = 1),
                        left_annotation=rowAnnotation(df=annotations, col=anno_colors, show_annotation_name = F))


pdf(paste0("outputs/figure_2/csx_closest_normal_heatmap_all_samples.pdf"), width = 5, height = 6)
print(plt)
dev.off()

# subset to stable samples
bulk_samples<-sample_info$ID[sample_info$is_bulk==1]
bulk_sample_info<-sample_info[sample_info$ID%in%bulk_samples,]
genotype_groups<-genotyping_info[,c("Sample.ID","Genomic.Bin")]
bulk_sample_info$genotype_groups<-genotype_groups$Genomic.Bin[match(bulk_sample_info$patient,genotype_groups$Sample.ID)]


cibersort_closest_normal_sub<-cibersort_closest_normal[cibersort_closest_normal$patient%in%bulk_sample_info$patient[bulk_sample_info$genotype_groups=="Stable"],]

heatmap_input<-dcast(cibersort_closest_normal_sub, patient_type ~ timepoint, value.var = "closest_norm")
heatmap_input<-column_to_rownames(heatmap_input, "patient_type")
heatmap_input<-heatmap_input[order(heatmap_input$REL),]
heatmap_input<-heatmap_input[order(heatmap_input$DX),]

annotations<-data.frame(type=strsplit(rownames(heatmap_input),"_") %>% sapply('[',2))
annotations$type<-tolower(annotations$type)

heatmap_input<-heatmap_input[order(annotations$type),]
annotations<-annotations[order(annotations$type),,drop=F]

celltype_pal<-c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB")
names(celltype_pal)<-c("B/Plasma-like","CLP-like","CMP/LMPP-like","Erythroid-like","GMP-like","HSC-like","Mono/DC/Baso-like","T/NK-like")
anno_colors<-list(type=c(lsc="yellow",bulk="red",blast="darkblue"))

plt<-ComplexHeatmap::Heatmap(heatmap_input, 
                        col=celltype_pal,
                        show_row_names=T, 
                        show_column_names=T, 
                        cluster_rows = T,
                        cluster_columns = F,
                        heatmap_legend_param = list(title = "Number Overlaps"),
                        column_title="Projection Closest Normals",
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        rect_gp = gpar(col = "white", lwd = 1),
                        left_annotation=rowAnnotation(df=annotations, col=anno_colors, show_annotation_name = F))

pdf(paste0("outputs/figure_2/csx_closest_normal_heatmap_stable_samples.pdf"), width = 5, height = 6)
print(plt)
dev.off()
