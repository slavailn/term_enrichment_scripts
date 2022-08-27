library(pathfindR)
library(reactome.db)
library(KEGGREST)
library(ggplot2)
library(gage)
library(gageData)
library(annotate)
library(org.Hs.eg.db)
## This worflow deals with clustering and visualization of term enrichment.
## This example is based on term enrichment results produced by ClueGO.
## Cluster and visualize pathway output from ClueGO using pathfindR package.
## Modify this script as necessary for diferrent inputs

setwd("/path/to/working_directory")
list.files()

clueGO <- read.csv("docs/ClueGO_results.csv", header = T)
# Select Reactome pathways
react_df <- clueGO[grep("REACTOME_Pathways", clueGO$Ontology.Source),]
head(react_df)[,1:4]
names(react_df)
react_df$GOID <- gsub(':', '-', react_df$GOID)

## Get total gene counts for each of REACTOME pathways
xx <- as.list(reactomePATHID2EXTID)

## Select human pathways and get number of genes per pathway
pathways <- lapply(xx, length)
pathways <- pathways[grep("HSA", names(pathways))]
pathways <- data.frame(REACT_ID=names(unlist(pathways)),
                       NUM_GENES=unlist(pathways))

head(pathways)
dim(pathways)

react_df <- merge(react_df, pathways, by.x="GOID", by.y="REACT_ID")
head(react_df)[1:2,]
sum(duplicated(react_df))

## Calculate number fold enrichment as follows: 
## 1. Calculate the proportion of DE genes relative to total genes per pathway -
## Prop_DE = Num_DE_genes/All_genes_pathway
## 2. Calculate the proportion of all genes in the pathway relative to all genes -
## Prop_ALL = All_genes_pathway/All_genes
Prop_DE <- (react_df$Nr..Genes/react_df$NUM_GENES) * 100
Prop_ALL <- (react_df$NUM_GENES/25000) * 100
react_df$FOLD_ENRICHMENT <- Prop_DE/Prop_ALL
head(react_df)
names(react_df)

## Reformat this data frame for pathfindR
path_df <- data.frame(ID=react_df$GOID, Term_Description=react_df$GOTerm,
                      Fold_Enrichment=react_df$FOLD_ENRICHMENT,
                      occurence=react_df$Nr..Genes,
                      lowest_p=react_df$Term.PValue.Corrected.with.Benjamini.Hochberg,
                      highest_p=react_df$Term.PValue.Corrected.with.Benjamini.Hochberg,
                      Down_regulated=react_df$Genes.Cluster..1,
                      Up_regulated=react_df$Genes.Cluster..2)
head(path_df)

## Remove brackets in 2 last columns
path_df$Down_regulated <- gsub('\\[|\\]', '', path_df$Down_regulated)
path_df$Up_regulated <- gsub('\\[|\\]', '', path_df$Up_regulated)

## Cluster enriched terms
path_df <- path_df[order(path_df$highest_p),] 
path_clustered <- cluster_enriched_terms(head(path_df, n=100), 
                                         use_description = T)
write.csv(path_clustered,file="reactome_pathway_clustering.csv")

## Enrichment bubble plot without clustering
enrichment_chart(result_df = path_df, 
                 top_terms = 20) + xlim(0, 500)
ggsave("REACTOME_DOTPLOT_TOP20.tiff", device="tiff", units="in", 
       width=8, height=6)

## Plot bubble chart and split by clusters
enrichment_chart(path_clustered, 
                 plot_by_cluster = TRUE) +
  xlim(0, 500)
ggsave("REACTOME_DOTPLOT_CLUSTERS_TOP100.tiff", device="tiff", units="in", 
       width=8, height=6)

## Plot term-gene map
term_gene_graph(result_df = path_df, use_description = TRUE)
ggsave("REACTOME_TERM_GENE_TOP100.tiff", device="tiff", units="in", 
       width=7, height=7)

## UpSet plot
tiff("REACTOME_UPSET_PLOT_TOP100.tiff", units="in", width=7, height=7,
     res=200)
UpSet_plot(result_df = path_df, num_terms = 15, use_description = T)
dev.off()

# ----------------------------------------------------------------- #
## Run the same visualization workflow for KEGG
## pathway analysis
data(kegg.sets.hs)
kegg.sets.hs[[1]]

# Translate entrez ids to gene symbols in the pathway list
kegg_hs <- 
  lapply(kegg.sets.hs, function(x) as.character(
    getSYMBOL(x, data='org.Hs.eg.db')))

pathways <- lapply(kegg_hs, length)
pathways <- data.frame(KEGG_ID=names(unlist(pathways)),
                       NUM_GENES=unlist(pathways))
head(pathways)[1:2,]
pathways$KEGG_ID <- gsub('^hsa| .+', '', pathways$KEGG_ID)
pathways$KEGG_ID <- paste("KEGG:", pathways$KEGG_ID, sep="")
rownames(pathways) <- pathways$KEGG_ID
head(pathways)

kegg_df <- clueGO[grep("KEGG", 
                       clueGO$Ontology.Source),]
clueGO$Ontology.Source[grep("KEGG", clueGO$Ontology.Source)]
head(kegg_df)[,1:6]

kegg_df <- merge(kegg_df, pathways, by.x="GOID", by.y="KEGG_ID")
head(kegg_df)[1:2,]
sum(duplicated(kegg_df))

## Calculate number fold enrichment as follows: 
## 1. Calculate the proportion of DE genes relative to total genes per pathway -
## Prop_DE = Num_DE_genes/All_genes_pathway
## 2. Calculate the proportion of all genes in the pathway relative to all genes -
## Prop_ALL = All_genes_pathway/All_genes
Prop_DE <- (kegg_df$Nr..Genes/kegg_df$NUM_GENES) * 100
Prop_ALL <- (kegg_df$NUM_GENES/25000) * 100
kegg_df$FOLD_ENRICHMENT <- Prop_DE/Prop_ALL
head(kegg_df)
names(kegg_df)

## Reformat this data frame for pathfindR
path_df <- data.frame(ID=kegg_df$GOID, Term_Description=kegg_df$GOTerm,
                      Fold_Enrichment=kegg_df$FOLD_ENRICHMENT,
                      occurence=kegg_df$Nr..Genes,
                      lowest_p=kegg_df$Term.PValue.Corrected.with.Benjamini.Hochberg,
                      highest_p=kegg_df$Term.PValue.Corrected.with.Benjamini.Hochberg,
                      Down_regulated=kegg_df$Genes.Cluster..1,
                      Up_regulated=kegg_df$Genes.Cluster..2)
head(path_df)

## Remove brackets in 2 last columns
path_df$Down_regulated <- gsub('\\[|\\]', '', path_df$Down_regulated)
path_df$Up_regulated <- gsub('\\[|\\]', '', path_df$Up_regulated)

## Cluster enriched terms
path_clustered <- cluster_enriched_terms(path_df, 
                                         use_description = T)
write.csv(path_clustered,file="kegg_pathway_clustering.csv")

## Enrichment bubble plot without clustering
enrichment_chart(result_df = path_df, 
                 top_terms = 20) + xlim(0, 50)
ggsave("KEGG_DOTPLOT_TOP20.tiff", device="tiff", units="in", 
       width=8, height=6)

## Plot bubble chart and split by clusters
enrichment_chart(path_clustered, plot_by_cluster = TRUE) +
  xlim(0, 50)
ggsave("KEGG_DOTPLOT_CLUSTERS.tiff", device="tiff", units="in", 
       width=8, height=6)

## Plot term-gene map
term_gene_graph(result_df = path_df, use_description = TRUE)
ggsave("KEGG_TERM_GENE.tiff", device="tiff", units="in", 
       width=7, height=7)

## UpSet plot
tiff("KEGG_UPSET_PLOT.tiff", units="in", width=7, height=7,
     res=200)
UpSet_plot(result_df = path_df, num_terms = 20, 
           use_description = T)
dev.off()
