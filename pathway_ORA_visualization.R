library(clusterProfiler)
library(pcaExplorer)
library(DOSE)
library(enrichplot)
library(GeneTonic)
library(DESeq2)
library(macrophage)
library("AnnotationDbi")
library(topGO)
library(GeneSetCluster)
library(readxl)
library(limma)
library(ggplot2)

## Clustering, visualization of enriched terms with pcaExplorer, GeneTonic, and ClusterProfiler

setwd("path/to/working_directory")
list.files()

data("geneList")

# ClusterProfiler vignette
de <- names(geneList)[abs(geneList) > 2]
edo <- enrichDGN(de)

edo2 <- gseDO(geneList)
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

## GeneTonic vignette
data("gse", package = "macrophage")
dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
dds_macrophage

keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
dds_macrophage

dds_macrophage <- DESeq(dds_macrophage)
res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
                                        contrast = c("condition", "IFNg", "naive"),
                                        lfcThreshold = 1, alpha = 0.05)
res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL

de_symbols_IFNg_vs_naive <- deseqresult2df(res_macrophage_IFNg_vs_naive, 
                                           FDR = 0.05)$SYMBOL

bg_ids <- rowData(dds_macrophage)$SYMBOL[rowSums(counts(dds_macrophage)) > 0]

topgoDE_macrophage_IFNg_vs_naive <-
  pcaExplorer::topGOtable(de_symbols_IFNg_vs_naive,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "symbol",
                          topTablerows = 500)

data("res_enrich_macrophage")
head(topgoDE_macrophage_IFNg_vs_naive, 2)

# Convert enrichment results to GeneTonic compatible format
res_enrich_macrophage <- 
  shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
colnames(res_enrich_macrophage)

# Get annotation data object for ID conversion
anno_df <- data.frame(
  gene_id = rownames(dds_macrophage),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds_macrophage), 
                     column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds_macrophage)
)
head(anno_df)

gs_alluvial(res_enrich = res_enrich_macrophage,
            res_de = res_macrophage_IFNg_vs_naive,
            annotation_obj = anno_df,
            n_gs = 4)
head(res_enrich_macrophage)
head(res_macrophage_IFNg_vs_naive)

# --------------------------------------------------------- #
# Analyze pathway enrichment data with GeneSetCluster
IPA.files <- c(system.file("extdata", "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
               system.file("extdata", "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
               system.file("extdata", "MM10.IPA.KO.uGvsMac.Functional_annotations.xls", package = "GeneSetCluster"),
               system.file("extdata", "MM10.IPA.WT.uGvsMac.Functional_annotations.xls", package = "GeneSetCluster"))
canonical.files <- IPA.files[grep("Canonical", IPA.files)]
canonical.files

MM10.IPA.KO.uGvsMac.Canonical <- read_excel(path = system.file("extdata", "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
                                            skip=1, sheet = 1)
MM10.IPA.WT.uGvsMac.Canonical <- read_excel(path = system.file("extdata", "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
                                            skip=1, sheet = 1)

MM10.IPA.KO.uGvsMac.Canonical <- 
  as.data.frame(MM10.IPA.KO.uGvsMac.Canonical)
MM10.IPA.WT.uGvsMac.Canonical <- 
  as.data.frame(MM10.IPA.WT.uGvsMac.Canonical)

head(MM10.IPA.KO.uGvsMac.Canonical)
head(MM10.IPA.WT.uGvsMac.Canonical)

# Calculating the number of molecules:
# we can see that the string is comma seperated for these molecules:
MM10.IPA.KO.uGvsMac.Canonical$MoleculesCount <- NA
head(MM10.IPA.KO.uGvsMac.Canonical)

for(can.i in 1:nrow(MM10.IPA.KO.uGvsMac.Canonical))
{
  mol.i <- 
    as.vector(strsplit2(as.character(MM10.IPA.KO.uGvsMac.Canonical[can.i,"Molecules"]), 
                               split=","))
  MM10.IPA.KO.uGvsMac.Canonical[can.i,"MoleculesCount"]<- length(mol.i)
}

head(MM10.IPA.KO.uGvsMac.Canonical)

MM10.IPA.KO.uGvsMac.Canonical.filtered <- MM10.IPA.KO.uGvsMac.Canonical[MM10.IPA.KO.uGvsMac.Canonical$`-log(p-value)` > 1.31 & 
                                                                          MM10.IPA.KO.uGvsMac.Canonical$MoleculesCount > 5,]

nrow(MM10.IPA.KO.uGvsMac.Canonical.filtered)
#We can see that we have 53 Gene-Sets which are significant according to our definition.

#Repeat for WT 
MM10.IPA.WT.uGvsMac.Canonical$MoleculesCount <- NA
for(can.i in 1:nrow(MM10.IPA.WT.uGvsMac.Canonical))
{
  mol.i <- as.vector(strsplit2(as.character(MM10.IPA.WT.uGvsMac.Canonical[can.i,"Molecules"]), split=","))
  MM10.IPA.WT.uGvsMac.Canonical[can.i,"MoleculesCount"]<- length(mol.i)
}

MM10.IPA.WT.uGvsMac.Canonical.filtered <- MM10.IPA.WT.uGvsMac.Canonical[MM10.IPA.WT.uGvsMac.Canonical$`-log(p-value)` > 1.31 & 
                                                         MM10.IPA.WT.uGvsMac.Canonical$MoleculesCount > 5,]

nrow(MM10.IPA.KO.uGvsMac.Canonical.filtered)
nrow(MM10.IPA.WT.uGvsMac.Canonical.filtered)

# Now we combine
# -> Pathways are just concatenated
# -> Molecules (aka the genes) are just concatenated,
# groups is a string that is the length of the combined pathways with the repeating info.
# -> Source is how the data was generated (for meta data reasons, not necessary to add)
# -> Type is what kind of data is it (for meta data reasons, not necessary to add)
# -> Structure is how the genes are presented, only important if you want to 
# combine gene sets, the genes have to match, 
# so the program wants to know its speaking the same language
# organism, same as the structure, only important for combining gene sets (optional)
# sep, how the genes in the molecules group are seperated. 
# -> Important for readign the individual genes.

IPA.KOvsWT.PathwayObject <- ObjectCreator(Pathways = c(MM10.IPA.KO.uGvsMac.Canonical.filtered$`Ingenuity Canonical Pathways`,
                                                       MM10.IPA.WT.uGvsMac.Canonical.filtered$`Ingenuity Canonical Pathways`), 
                                          Molecules = c(MM10.IPA.KO.uGvsMac.Canonical.filtered$Molecules,
                                                        MM10.IPA.WT.uGvsMac.Canonical.filtered$Molecules),
                                          Groups = c(rep("KO", times = nrow(MM10.IPA.KO.uGvsMac.Canonical.filtered)),
                                                     rep("WT", times = nrow(MM10.IPA.WT.uGvsMac.Canonical.filtered))),
                                          Source = "IPA",
                                          Type = "Canonical_Pathways",#Optional
                                          structure = "SYMBOL",
                                          organism ="org.Mm.eg.db",
                                          sep = ",")
IPA.KOvsWT.PathwayObject

# Create background object
Great.files <- c(system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed.tsv", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed.tsv", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"))
Great.files.bckgrnd <- Great.files[grepl("BCKGRND", Great.files)]

Great.bckgnrd.Object1 <- LoadGeneSets(file_location = Great.files.bckgrnd, 
                                      groupnames= c("KO", "WT"),
                                      P.cutoff = 0.05, 
                                      Mol.cutoff = 5,
                                      Source = "Great",
                                      Great.Background = T,#specify the background, as great has a different output if run with or without background
                                      type = "Canonical_Pathways",
                                      topranks = 20,#Great gives soo much output, recommended is adding a topranks filter for first 20
                                      structure = "SYMBOL",
                                      Organism = "org.Mm.eg.db",
                                      seperator = ",")

## Select pathway types of interest
man.Great.Object1 <- ManageGeneSets(Object = Great.bckgnrd.Object1, 
                                    keep.type = c("Disease Ontology","GO Biological Process" ), 
                                    exclude.type = "")


ShowExperimentdata(Object = man.Great.Object1 )
ShowMeta(Object = man.Great.Object1 )

## Combine and cluster gene sets
man.Great.Object2 <- CombineGeneSets(Object = man.Great.Object1)

## Determine optimal number of clusters
OptimalGeneSets(object = man.Great.Object2, method = "silhouette", max_cluster= 24, 
                cluster_method = "kmeans", main= "Kmeans for 24 clusters")


## Cluster gene sets
man.Great.Object3 <- ClusterGeneSets(Object = man.Great.Object2, 
                                     clusters = 5, 
                                     method = "kmeans")

## Plot gene sets
## K means without scaling
PlotGeneSets(Object = man.Great.Object3, fontsize = 3,
             legend = T,
             annotation.mol = F,
             main="Great_Background clustered with Kmeans without scaling \n Disease Ontology and GO Biological Process")

## K-means without scaling
PlotGeneSets(Object =man.Great.Object3, fontsize = 3,
             legend = T,
             annotation.mol=F,
             RR.max = 60,
             main="Great_Background clustered with Kmeans \n Disease Ontology and GO Biological Process")
