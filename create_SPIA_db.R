library(SPIA)
library(pathview)


setwd(<wd>)

# Download KEGG pathway xml file with KEGGREST
# This example uses Rattus norvegicus
listDatabases()
org <- keggList("organism")
head(org)
kegg_rno <- keggList("pathway", "rno")
kegg_rno <- names(kegg_rno)
kegg_rno <- gsub("rno", "", kegg_rno)

download.kegg(pathway.id = kegg_rno, species = "rno", 
              kegg.dir = ".",
              file.type=c("xml"))

# Make SPIA data for rat
makeSPIAdata(kgml.path="../kegg_xml",organism="rno", out.path=".")



