# Armon Azizi (aazizi@stanford.edu)
# 
# Nuno/Azizi et al.
# Figure 4
# comparison of projection fraction celltypes and bulk cibersort deconvolution fractions
#

library(ArchR)
library(openxlsx)
library(Seurat)

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")


cibersort_fractions<-as.data.frame(fread("~/Bioinformatics/AML_relapse_project/cibersort/deconvolution_output/AML_vs_normal_heme_signature/CIBERSORTx_Results.txt"))

scATAC_ref<-readRDS("inputs/scATAC/sample_reference_hg19.rds")

sample_info<-openxlsx::read.xlsx("inputs/bulk_sample_info/sample_info.xlsx")

# compute single cell projection fractions
AML_projection_fractions<-data.frame()
for(p in unique(scATAC_ref$patient)){
  assignments<-readRDS(paste0("outputs/scATAC/AML_scATAC_projection/",p,"_LSI_projection.rds"))
  assignments<-assignments[assignments$AML_classification!="reference",]
  assignment_counts<-as.data.frame(table(assignments[,c("Sample","AML_classification")]))
  for(s in unique(assignment_counts$Sample)){assignment_counts$Freq[assignment_counts$Sample==s]<-assignment_counts$Freq[assignment_counts$Sample==s]/sum(assignment_counts$Freq[assignment_counts$Sample==s])}
  AML_projection_fractions<-rbind(AML_projection_fractions, data.frame(patient=assignment_counts$Sample,
                                                                       classification=assignment_counts$AML_classification,
                                                                       frequency=assignment_counts$Freq))
}


cibersort_fractions$`B/Plasma`<-rowSums(cibersort_fractions[,c("Bcell"),drop=F])
cibersort_fractions$CLP<-rowSums(cibersort_fractions[,c("CLP"),drop=F])
cibersort_fractions$`CMP/LMPP`<-rowSums(cibersort_fractions[,c("LMPP", "CMP")])
cibersort_fractions$Erythroid<-rowSums(cibersort_fractions[,c("Ery","MEP")])
cibersort_fractions$GMP<-rowSums(cibersort_fractions[,c("GMP"),drop=F])
cibersort_fractions$HSC<-rowSums(cibersort_fractions[,c("HSC"),drop=F])
cibersort_fractions$`Mono`<-rowSums(cibersort_fractions[,c("Mono"),drop=F])
cibersort_fractions$`T/NK`<-rowSums(cibersort_fractions[,c("CD4Tcell","CD8Tcell","NKcell")])
cibersort_fractions$sample<-sample_info$sample[match(sample_info$name, cibersort_fractions$Mixture)]
cibersort_fractions$is_bulk<-sample_info$is_bulk[match(sample_info$name, cibersort_fractions$Mixture)]

for(i in 1:nrow(AML_projection_fractions)){
  p<-AML_projection_fractions[i,"patient"]
  classification<-as.character(AML_projection_fractions[i,"classification"])
  cibersort_freq<-mean(cibersort_fractions[cibersort_fractions$sample==p&cibersort_fractions$is_bulk==1,sapply(strsplit(classification,"-"),'[',1)])
  AML_projection_fractions[i,"cibersort_freq"]<-cibersort_freq
}

AML_projection_fractions<-AML_projection_fractions[complete.cases(AML_projection_fractions),]

celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono-like","DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))

plt<-ggplot(AML_projection_fractions, aes(x=log10(frequency), y=log10(cibersort_freq), color=classification)) +
  geom_point(size=2) +
  theme_classic() +
  scale_color_manual(values=celltype_pal, name = "Cell Type") +
  scale_x_continuous(limits = c(-4,1)) +
  scale_y_continuous(limits = c(-4,1)) +
  xlab("Log10(Single Cell Projection Frequency)") +
  ylab("Log10(CSX Fraction)") +
  stat_cor(method = "pearson", mapping = aes(x=frequency, y=cibersort_freq), inherit.aes = F) +
  geom_smooth(method = "lm", se = FALSE, color="red")

pdf(paste0("outputs/scATAC_plots/cibersort_fractions_vs_projection_classification.pdf"), width = 5, height = 3.5, useDingbats = FALSE)
print(plt)
dev.off()

# ggplot(AML_projection_fractions, aes(x=log10(frequency), y=log10(cibersort_freq), color=patient)) +
#   geom_point() +
#   theme_classic() +
#   scale_x_continuous(limits = c(-4,1)) +
#   scale_y_continuous(limits = c(-4,1)) +
#   stat_cor(method = "pearson", mapping = aes(x=frequency, y=cibersort_freq), inherit.aes = F) +
#   stat_cor(method = "pearson") +
#   geom_smooth(method = "lm", se = FALSE, color="red")
# 
# ggplot(AML_projection_fractions, aes(x=frequency, y=cibersort_freq, color=classification)) +
#   geom_point() +
#   theme_classic() +
#   stat_cor(method = "pearson", mapping = aes(x=frequency, y=cibersort_freq), inherit.aes = F) +
#   geom_smooth(method = "lm", se = FALSE, color="red")
# 
# ggplot(AML_projection_fractions, aes(x=frequency, y=cibersort_freq, color=patient)) +
#   geom_point() +
#   theme_classic() +
#   stat_cor(method = "pearson", mapping = aes(x=frequency, y=cibersort_freq), inherit.aes = F) +
#   geom_smooth(method = "lm", se = FALSE, color="red")


## bulk classifications and comparison of closest normal celltype
bulk_classifications<-readRDS("healthy_hematopoiesis_plus_AML/AML_cell_classifications/bulk_KNN_classifications.rds")

classification_celltypes<-c("HSC","CMP/LMPP","CLP","GMP","B/Plasma","Erythroid","Mono/DC/Baso", "T/NK")

cibersort_closest_normal<-data.frame(mixture=cibersort_fractions$Mixture,
                                     closest_norm=classification_celltypes[max.col(cibersort_fractions[,classification_celltypes],ties.method="first")])

# get closest norm for bulk projections
bulk_classifications<-bulk_classifications[bulk_classifications$type!="scATAC",]
aml_samples<-unique(bulk_classifications$sample)
projection_closest_normal<-data.frame(sample=aml_samples,
                                      closest_norm=sapply(aml_samples, function(s) names(sort(table(bulk_classifications$AML_classification[bulk_classifications$sample==s]),decreasing=TRUE)[1])))

closest_norms<-data.frame(sample=aml_samples,
                          projection_closest_norm=projection_closest_normal$closest_norm,
                          cibersort_closest_norm=cibersort_closest_normal$closest_norm[match(aml_samples, cibersort_closest_normal$mixture)])
closest_norms$projection_closest_norm<-strsplit(closest_norms$projection_closest_norm,"-") %>% sapply('[',1)

overlaps<-as.data.frame(table(closest_norms[,c("projection_closest_norm","cibersort_closest_norm")]))
overlaps<-dcast(overlaps, projection_closest_norm ~ cibersort_closest_norm)
overlaps<-column_to_rownames(overlaps, "projection_closest_norm")

overlaps<-overlaps[classification_celltypes[classification_celltypes%in%rownames(overlaps)],classification_celltypes[classification_celltypes%in%colnames(overlaps)]]

ComplexHeatmap::Heatmap(overlaps, 
                        #col=viridis::viridis(50),
                        col=colorRamp2(seq(0,50,length.out = 10), viridis::viridis(10)),
                        show_row_names=T, 
                        show_column_names=T, 
                        cluster_rows = F,
                        cluster_columns = F,
                        heatmap_legend_param = list(title = "Number Overlaps"),
                        column_title="Cibersort Classification", 
                        column_title_side = "bottom",
                        row_names_side = "left",
                        row_title = "Projection Classification",
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 10))



## summarize projected closest celltypes
bulk_classifications<-readRDS("healthy_hematopoiesis_plus_AML/AML_cell_classifications/bulk_KNN_classifications.rds")

classification_celltypes<-c("HSC","CMP/LMPP","CLP","GMP","B/Plasma","Erythroid","Mono/DC/Baso", "T/NK")

# get closest norm for bulk projections
bulk_classifications<-bulk_classifications[bulk_classifications$type!="scATAC",]
bulk_classifications$sample<-substr(bulk_classifications$sample,1,nchar(bulk_classifications$sample)-2)
aml_samples<-unique(bulk_classifications$sample)
projection_closest_normal<-data.frame(sample=aml_samples,
                                      closest_norm=sapply(aml_samples, function(s) names(sort(table(bulk_classifications$AML_classification[bulk_classifications$sample==s]),decreasing=TRUE)[1])))
projection_closest_normal$patient<-strsplit(rownames(projection_closest_normal),"_") %>% sapply('[',1)
projection_closest_normal$timepoint<-strsplit(rownames(projection_closest_normal),"_") %>% sapply('[',2)
projection_closest_normal$type<-strsplit(rownames(projection_closest_normal),"_") %>% sapply('[',3)
projection_closest_normal$patient_type<-paste(projection_closest_normal$patient,projection_closest_normal$type, sep="_")

heatmap_input<-dcast(projection_closest_normal, patient_type ~ timepoint, value.var = "closest_norm")
heatmap_input<-column_to_rownames(heatmap_input, "patient_type")
heatmap_input<-heatmap_input[order(heatmap_input$REL),]
heatmap_input<-heatmap_input[order(heatmap_input$DX),]

annotations<-data.frame(type=strsplit(rownames(heatmap_input),"_") %>% sapply('[',2))
annotations$type<-tolower(annotations$type)

heatmap_input<-heatmap_input[order(annotations$type),]
annotations<-annotations[order(annotations$type),,drop=F]

celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono/DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
anno_colors<-list(type=c(lsc="yellow",bulk="red",blast="darkblue"))

ComplexHeatmap::Heatmap(heatmap_input, 
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

# subset to stable samples
sample_info<-read.xlsx("../sample_info/sample_info.xlsx")
genotyping_info<-read.xlsx("../sample_info/genotyping_final.xlsx")
bulk_samples<-sample_info$ID[sample_info$is_bulk==1]
bulk_sample_info<-sample_info[sample_info$ID%in%bulk_samples,]
genotype_groups<-genotyping_info[,c("Sample.ID","Genomic.Bin")]
bulk_sample_info$genotype_groups<-genotype_groups$Genomic.Bin[match(bulk_sample_info$patient,genotype_groups$Sample.ID)]


projection_closest_normal_sub<-projection_closest_normal[projection_closest_normal$patient%in%bulk_sample_info$patient[bulk_sample_info$genotype_groups=="Stable"],]

heatmap_input<-dcast(projection_closest_normal_sub, patient_type ~ timepoint, value.var = "closest_norm")
heatmap_input<-column_to_rownames(heatmap_input, "patient_type")
heatmap_input<-heatmap_input[order(heatmap_input$REL),]
heatmap_input<-heatmap_input[order(heatmap_input$DX),]

annotations<-data.frame(type=strsplit(rownames(heatmap_input),"_") %>% sapply('[',2))
annotations$type<-tolower(annotations$type)

heatmap_input<-heatmap_input[order(annotations$type),]
annotations<-annotations[order(annotations$type),,drop=F]

celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono/DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
anno_colors<-list(type=c(lsc="yellow",bulk="red",blast="darkblue"))

ComplexHeatmap::Heatmap(heatmap_input, 
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



# cibersort closest normal summary
cibersort_fractions$sample<-substr(cibersort_fractions$Mixture,1,nchar(cibersort_fractions$Mixture)-2)

cibersort_fractions<-aggregate(cibersort_fractions[,classification_celltypes], by=list(cibersort_fractions$sample), FUN=mean)

cibersort_closest_normal<-data.frame(mixture=cibersort_fractions$Group.1,
                                     closest_norm=classification_celltypes[max.col(cibersort_fractions[,classification_celltypes],ties.method="first")])


cibersort_closest_normal$patient<-strsplit(cibersort_closest_normal$mixture,"_") %>% sapply('[',1)
cibersort_closest_normal$timepoint<-strsplit(cibersort_closest_normal$mixture,"_") %>% sapply('[',2)
cibersort_closest_normal$type<-strsplit(cibersort_closest_normal$mixture,"_") %>% sapply('[',3)
cibersort_closest_normal$patient_type<-paste(cibersort_closest_normal$patient,projection_closest_normal$type, sep="_")
cibersort_closest_normal$closest_norm<-paste0(cibersort_closest_normal$closest_norm, "-like")

heatmap_input<-dcast(cibersort_closest_normal, patient_type ~ timepoint, value.var = "closest_norm")
heatmap_input<-column_to_rownames(heatmap_input, "patient_type")
heatmap_input<-heatmap_input[order(heatmap_input$REL),]
heatmap_input<-heatmap_input[order(heatmap_input$DX),]


annotations<-data.frame(type=strsplit(rownames(heatmap_input),"_") %>% sapply('[',2))
annotations$type<-tolower(annotations$type)

heatmap_input<-heatmap_input[order(annotations$type),]
annotations<-annotations[order(annotations$type),,drop=F]

celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono/DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
celltype_pal<-celltype_pal[names(celltype_pal)%in%unique(unlist(heatmap_input))]
anno_colors<-list(type=c(lsc="yellow",bulk="red",blast="darkblue"))

ComplexHeatmap::Heatmap(heatmap_input, 
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


# subset to stable samples
sample_info<-read.xlsx("../sample_info/sample_info.xlsx")
genotyping_info<-read.xlsx("../sample_info/genotyping_final.xlsx")
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

celltype_pal<-ArchR::paletteDiscrete(c("HSC-like","CMP/LMPP-like","Mono/DC/Baso-like","GMP-like","CLP-like","B/Plasma-like","Erythroid-like","T/NK-like"))
anno_colors<-list(type=c(lsc="yellow",bulk="red",blast="darkblue"))

ComplexHeatmap::Heatmap(heatmap_input, 
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
