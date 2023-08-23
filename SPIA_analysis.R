library(SPIA)

# Set working directory
setwd(<wd>)

resFile <-  "path/to/DESeq2_results" # has to have entrez ids as annotations
compName <- "control_vs_treated" # Name of the comparison

res <- read.csv(resFile, header = T)
  
# Select differentially expressed genes (adjusted p-value < 0.05)
# This is a named vector, where names are entrez ids, 
# and values are log2 fold changes
ids <- res[res$padj < 0.05,]$entrezgene_id
logFC <- res[res$padj < 0.05,]$log2FoldChange
DE_genes <- logFC
names(DE_genes) <- ids
DE_genes <- DE_genes[!is.na(names(DE_genes))] # Remove NAs
DE_genes <- DE_genes[!duplicated(names(DE_genes))] # Remove duplicates
  
# Get all expressed genes entrez ids
ALL_genes <- res$entrezgene_id
ALL_genes <- ALL_genes[!is.na(ALL_genes)] # remove NAs
ALL_genes <- ALL_genes[!duplicated(ALL_genes)] # remove duplicates
res <- spia(de=DE_genes, all=ALL_genes, organism="rno", 
            nB=2000, plots=FALSE, beta=NULL, combine="fisher", 
            verbose=FALSE, data.dir = "spia_db/") # point to folder with SPIA database
write.csv(res, file=paste(compName, "_SPIA_result.csv", sep = ""))
tiff(file=paste(compName, "_SPIA_result.tiff", sep = ""), 
      width=500, height=500)
plotP(res, threshold = 0.05)
dev.off()
