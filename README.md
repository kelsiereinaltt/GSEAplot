# GSEAplot R package

Gene Set Enrichment Analysis (GSEA) is used to identify differentially expressed gene sets that are enriched for annotated biological functions. The existing GSEA R code is not in the form of a flexible package with analysis and plotting customization options, and the results produced are not generated in the form of R objects. Here we introduce the GSEAplot R package for saving relevant information from the analysis to the current R workspace, as well as introducing the ability to customize plots and databases. The GSEAplot package provides a novel utility that facilitates the implementation of  GSEA in genomics analysis pipelines. 

## Documentation

See the vignette and user manual:  
GSEAplot.vignette.pdf  
GSEAplot-manual.pdf

## Installation

devtools::install_github("kelsiereinaltt/GSEAplot")  
library(GSEAplot)

## Authors

Sarah Innis, Kelsie Reinaltt, Mete Civelek, Warren Anderson
