library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)
library(enrichplot)
library(GOSemSim)
library(ggtree)
library(ggupset)

setwd(<working_directory>)

res <- read.csv("HGHL_vs_CT_PSI/HGHL_vs_CT_PSI_results.csv", header = T)
head(res)[,1:3]

# select significant entries
res_sig <- res[res$padj < 0.05,]
res_sig <- res_sig$entrezgene_id
res_sig <- res_sig[!is.na(res_sig)]
res_sig

# Universe ids
universe_ids <- res$entrezgene_id[!is.na(res$entrezgene_id)]
universe_ids

# Run enrichment analysis
ego <- enrichGO(gene          = res_sig,
                universe      = as.character(universe_ids),
                OrgDb         = org.Rn.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write.csv(ego@result, file="HGHL_vs_CT_PSI_GO_clusterProfiler.csv")

dotplot(ego, showCategory=30) + ggtitle("")
ggsave("HGHL_vs_CT_PSI_GO_dotplot.pdf", units = "in", width = 6, 
       height = 12)
d <- godata('org.Rn.eg.db', ont="BP")
ego2 <- pairwise_termsim(ego, method = "Wang", semData = d)
enrichplot::treeplot(ego2, method = "average")
ggsave("HGHL_vs_CT_PSI_GO_treeplot.pdf", units = "in", width = 12, 
       height = 12)

# Build GO graph for GO enrichment results
goplot(ego, showCategory = sum(ego@result$p.adjust < 0.05))
ggsave("HGHL_vs_CT_PSI_GO_graph.pdf", units = "in", width = 15, 
       height = 15)

# Build gene to gene-set network
# Get fold changes
geneList <- res$log2FoldChange
names(geneList) <- res$entrezgene_id
head(geneList)
geneList <- geneList[!is.na(names(geneList))]
egox <- setReadable(ego, 'org.Rn.eg.db', 'ENTREZID')
p1 <- cnetplot(egox, color.params = list(foldChange = geneList),
               showCategory = 10)
p1
ggsave("HGHL_vs_CT_PSI_gene_set_graph.pdf", units = "in", width = 20, 
       height = 20)

# Show enriched categories and genes on the heatmap
p2 <- heatplot(egox, foldChange=geneList, showCategory=10)
p2
ggsave("HGHL_vs_CT_GO_heatmap.pdf", units = "in", width = 20, 
       height = 5)

# Build enrichment map
ego <- pairwise_termsim(ego)
p3 <- emapplot(ego, cex_category=1.5, layout="kk", 
               showCategory=sum(ego@result$p.adjust < 0.05))
p3
ggsave("HGHL_vs_CT_PSI_GO_emaplot.pdf", units = "in", width = 20, 
       height = 20)


##########################################################################
## Run KEGG enrichment analysis
res <- read.csv("HGHL_vs_CT_PSI/HGHL_vs_CT_PSI_results.csv", header = T)

# select significant entries
res_sig <- res[res$padj < 0.05,]
res_sig <- res_sig$entrezgene_id
res_sig <- res_sig[!is.na(res_sig)]
res_sig

# Create a vector of universe ids
universe_ids <- res$entrezgene_id[!is.na(res$entrezgene_id)]
universe_ids

# Kegg enrichment
kk <- enrichKEGG(gene         = as.character(res_sig),
                 universe = as.character(universe_ids),
                 organism     = 'rno',
                 pvalueCutoff = 0.05,
                 )
head(kk)
kk@result
write.csv(kk@result, file="HGHL_vs_CT_PSI_KEGG_clusterProfiler.csv")

# Dotplot for KEGG enrichment
dotplot(kk, showCategory=11) + ggtitle("dotplot for KEGG enrichment")
ggsave("HGHL_vs_CT_PSI_KEGG_dotplot.pdf", units = "in", width = 6, 
       height = 5)

# Build Gene-Concept network
kk_entrez <- setReadable(kk, 'org.Rn.eg.db', 'ENTREZID')
geneList <- res$log2FoldChange
names(geneList) <- res$entrezgene_id
head(geneList)
geneList <- geneList[!is.na(names(geneList))]
p1 <- cnetplot(kk_entrez, foldChange=geneList,
               showCategory = 11)
p1
ggsave("HGHL_vs_CT_PSI_KEGG_cnet_plot.pdf", units = "in", width = 15, 
       height = 15)

# build a heatplot
p2 <- heatplot(kk_entrez, foldChange=geneList, showCategory=11)
p2
ggsave("HGHL_vs_CT_PSI_KEGG_heatplot.pdf", units = "in", width = 30, 
       height = 6)

## ------------------------------------------------------------ ##
## Run wikipathways enrichment analysis
## Only one significantly enriched pathway
wp <- enrichWP(as.character(res_sig), universe = as.character(universe_ids),
                pvalueCutoff = 0.05,
                organism = "Rattus norvegicus")

wp@result

write.csv(wp@result, file="HGHL_vs_CT_PSI_WP_clusterProfiler.csv")

# Dotplot for KEGG enrichment
dotplot(wp, showCategory=1) + ggtitle("dotplot for WikiPathways enrichment")
ggsave("HGHL_vs_CT_PSI_WP_dotplot.pdf", units = "in", width = 6, 
       height = 5)

# Build Gene-Concept network
wp_entrez <- setReadable(wp, 'org.Rn.eg.db', 'ENTREZID')
geneList <- res$log2FoldChange
names(geneList) <- res$entrezgene_id
head(geneList)
geneList <- geneList[!is.na(names(geneList))]
p1 <- cnetplot(wp_entrez, foldChange=geneList,
               showCategory = 1)
p1
ggsave("HGHL_vs_CT_PSI_WP_cnet_plot.pdf", units = "in", width = 10, 
       height = 10)


# build a heatplot
p2 <- heatplot(wp_entrez, foldChange=geneList, showCategory=1)
p2
ggsave("HGHL_vs_CT_PSI_WP_heatplot.pdf", units = "in", width = 6, 
       height = 3)
