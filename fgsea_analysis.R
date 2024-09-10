## Find enriched GO categories using fgsea
## In this example we are using fgsea to calculate enrichment
## relative to GO slim
## The genes were ranked based on p-values
library(fgsea)
library(data.table)
library(ggplot2)
library(biomaRt)
library(GO.db)

# Download GO slim categories
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")
go.slim <- getBM(attributes = "goslim_goa_accession", mart = mart)[,1]
results <- AnnotationDbi::select(org.Hs.eg.db, keys = go.slim, 
                                 columns = c("SYMBOL"), 
                                 keytype = "GOALL")
head(results)

# Select GO categories biological process
# Create a list of vectors where names are 
results_bp <- subset(results, ONTOLOGYALL == "BP")
head(results_bp)
go_bp <- split(results_bp$SYMBOL, results_bp$GOALL)
head(go_bp)

# Convert entrez gene ids in meta-analysis results to gene symbols
getwd()
list.files()
wfish <- read.csv("meta_analysis_results/meta_analysis_wfisher.csv",
                  header = T)
head(wfish)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
attributes <- c("entrezgene_id", "hgnc_symbol")
ids <- wfish$entrezID
annot <- getBM(attributes = attributes, values = as.character(ids), 
               mart = ensembl)
annot <- annot[!duplicated(annot$entrezgene_id),]
annot <- annot[!is.na(annot$entrezgene_id),]
head(annot)

wfish_annot <- merge(wfish, annot, by.x = "entrezID", 
                     by.y = "entrezgene_id")
head(wfish_annot)
sum(duplicated(wfish_annot$entrezID))
sum(duplicated(wfish_annot$hgnc_symbol))
wfish_annot <- wfish_annot[!duplicated(wfish_annot$hgnc_symbol),]

# Rank meta-analysis result by p-value
wfish_annot <- wfish_annot[order(wfish_annot$pvalue),]
ranked_list <- wfish_annot$pvalue
names(ranked_list) <- wfish_annot$hgnc_symbol
head(ranked_list)

# Run GSEA analysis relative to GO slim "Biological process"
# Run fgsea
fgseaRes <- fgsea(pathways = go_bp, 
                  stats    = ranked_list,
                  minSize  = 15,
                  maxSize  = 500)
head(fgseaRes[order(pval), ], n=10)
plotEnrichment(go_bp[["GO:0007018"]],
               ranked_list) + labs(title="GO:0007018")

# Make a table plot for a number of selected pathways
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
plotGseaTable(go_bp[topPathways], ranked_list, fgseaRes, 
              gseaParam=0.5)
