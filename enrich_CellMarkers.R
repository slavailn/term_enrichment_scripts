library(org.Hs.eg.db)

# Calculate gene enrichment relative to cell types
setwd(<wd>)
list.files()

wfisher_res <- read.csv("meta_analysis_results/meta_analysis_wfisher.csv", 
                        header = T)
head(wfisher_res)
# Get a list of all up-regulated genes with p.value below 0.01
wfisher_res <- subset(wfisher_res, pvalue < 0.01 & effect_size_dir == '+')
dim(wfisher_res)

# Add hgnc symbols
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
xx <- unlist(xx)
annot <- data.frame(entrez_id = names(xx), hgnc_symbol = xx)
head(annot)
wfisher_annot <- merge(wfisher_res, annot, by.x = "entrezID", 
                       by.y = "entrez_id")
head(wfisher_annot)
sum(duplicated(wfisher_annot$entrezID))

write.csv(wfisher_annot, file="wfisher_up_p0.01.csv")

# Load the list of Up-regulated genes with p-value less than 0.01
# into WebCSEA and download the results https://bioinfo.uth.edu/webcsea/
