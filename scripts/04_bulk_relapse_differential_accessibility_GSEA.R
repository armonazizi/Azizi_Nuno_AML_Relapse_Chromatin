# Armon Azizi (aazizi@stanford.edu)
# Nuno/Azizi et al.
# Figure 2
# GSEA of differential gene accessibility in stable relapse samples

setwd("~/Bioinformatics/AML_relapse_project/analysis_clean/")
source("scripts/00_global_ATAC_functions.R")

pathways <- c(gmtPathways("inputs/GSEA/c2.all.v7.2.symbols.gmt"),
              gmtPathways("inputs/GSEA/h.all.v7.2.symbols.gmt"))

# Do gsea on all pathways for stable samples
ranked_genes<-readRDS("outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")
ranked_genes <- setNames(ranked_genes$log2FoldChange, rownames(ranked_genes))
ranked_genes<-ranked_genes[!is.na(ranked_genes)]

res_stable<-fgsea(pathways, ranked_genes)

# Do gsea on all pathways for nonstable samples
ranked_genes<-readRDS("outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_nonstable_cases.rds")
ranked_genes <- setNames(ranked_genes$log2FoldChange, rownames(ranked_genes))
ranked_genes<-ranked_genes[!is.na(ranked_genes)]
res_nonstable<-fgsea(pathways, ranked_genes)


# Combine all data and save to excel
res_combined<-data.frame(pathway=res_stable$pathway,
                         stable_pval=res_stable$pval,
                         stable_padj=res_stable$padj,
                         stable_log2err=res_stable$log2err,
                         stable_ES=res_stable$ES,
                         stable_NES=res_stable$NES,
                         stable_size=res_stable$size,
                         nonstable_pval=res_nonstable$pval[match(res_stable$pathway, res_nonstable$pathway)],
                         nonstable_padj=res_nonstable$padj[match(res_stable$pathway, res_nonstable$pathway)],
                         nonstable_log2err=res_nonstable$log2err[match(res_stable$pathway, res_nonstable$pathway)],
                         nonstable_ES=res_nonstable$ES[match(res_stable$pathway, res_nonstable$pathway)],
                         nonstable_NES=res_nonstable$NES[match(res_stable$pathway, res_nonstable$pathway)],
                         nonstable_size=res_nonstable$size[match(res_stable$pathway, res_nonstable$pathway)],
                         geneset_group=rep("msigdb_curated", nrow(res_stable)))
res_combined$geneset_group[res_combined$pathway%in%names(gmtPathways("inputs/GSEA/h.all.v7.2.symbols.gmt"))]<-"msigdb_hallmark"
res_combined$genes<-res_stable$leadingEdge

write.xlsx(res_combined, "outputs/bulk_analysis/bulk_stable_and_nonstable_gene_scores_gsea_res.xlsx")


# generate enrichment plots for select pathways in stable samples
ranked_genes<-readRDS("outputs/bulk_analysis/REL_vs_DX_gene_accessibility_deseq_res_stable_cases.rds")
ranked_genes <- setNames(ranked_genes$log2FoldChange, rownames(ranked_genes))
ranked_genes<-ranked_genes[!is.na(ranked_genes)]

pathways_of_interest<-c("HALLMARK_MYC_TARGETS_V1",
                        "HALLMARK_E2F_TARGETS")
for(p in pathways_of_interest){
  res_sub<-fgsea(pathways[p], ranked_genes)
  pval<-signif(res_sub$padj, digits=3)
  nes<-signif(res_sub$NES, digits=3)
  plt<-plotEnrichment(pathways[[p]],ranked_genes) +
    ggtitle(paste0(p, " Stable Samples")) +
    annotate(geom = 'text', label = paste0("p=",pval,"\nNES=",nes), x = -Inf, y = -Inf, hjust = 0, vjust = -3)
  
  print(plt)
  
  pdf(paste0("outputs/figure_2/GSEA/",p,"_stable_bulk.pdf"), width = 5, height = 4)
  print(plt)
  dev.off()
}
