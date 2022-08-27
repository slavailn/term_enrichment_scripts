library("rrvgo")
library("org.Hs.eg.db")

## Semantic space analysis of term enrichment results with rrvgo

setwd("path/to/working_directory")

## Vignette
go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))

simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
head(scores)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
head(reducedTerms)

## Plot heatmap based on similarity
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

## Plot scatter
scatterPlot(simMatrix, reducedTerms)

## Plot treemap
treemapPlot(reducedTerms)

## World cloud
wordcloudPlot(reducedTerms, min.freq=1, colors="black")

## Shiny app
rrvgo::shiny_rrvgo()

# --------------------------------------------------------- #
## Create scatter plots based on proteomic GOclue data
list.files()

goFile="Biol_Proc_Higher_WT.csv"
scatterName="Biol_Process_Higher_WT_scatter.tiff"
treemapName="Biol_Process_Higher_WT_treemap.tiff"
ontology="BP"

plotGOScatter <- function(goFile, scatterName, treemapName, ontology) 
{
  go_terms <- read.csv(goFile, header=T)
  go_terms <- go_terms[go_terms$Term.PValue.Corrected.with.Benjamini.Hochberg 
                       < 0.05,]
  
  simMatrix <- calculateSimMatrix(go_terms$GOID,
                                  orgdb="org.Mm.eg.db",
                                  ont=ontology,
                                  method="Rel")
  
  scores <- setNames(-log10(go_terms$Term.PValue.Corrected.with.Benjamini.Hochberg), 
                     go_terms$GOID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Mm.eg.db")
  p1 <- scatterPlot(simMatrix, reducedTerms) 
  ggsave(file=scatterName, plot=p1, device="tiff", width=6, height = 6,
         units = "in")
  
  tiff(file=treemapName, res=300, units = "in", width = 6, height = 6)
  treemapPlot(reducedTerms)
  dev.off()
  
  
}


plotGOScatter(goFile="Biol_Proc_Higher_WT.csv", 
              scatterName="Biol_Process_Higher_WT_scatter.tiff",
              treemapName="Biol_Process_Higher_WT_treemap.tiff",
              ontology="BP")

plotGOScatter(goFile="Cell_Comp_Higher_WT.csv", 
              scatterName="Cell_Comp_Higher_WT_scatter.tiff",
              treemapName="Cell_Comp_Higher_WT_treemap.tiff",
              ontology="CC")

plotGOScatter(goFile="Mol_Func_Higher_WT.csv", 
              scatterName="Mol_Func_Higher_WT_scatter.tiff",
              treemapName="Mol_Func_Higher_WT_treemap.tiff",
              ontology="MF")















