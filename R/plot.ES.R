#' GSEA Enrichment Plots
#'
#' The plot.ES function creates a pdf of enrichment plots obtained from the GSEA analysis given a list of ggplot objects and a character by which to name the pdf. The pdf is formatted to contain 4 graphs per page and save to the working directory.
#' @param list.of.plots A list of ggplot objects containing the necessary information for an enrichment plot using ggplot2.
#' @param plotname A character by which to name and save the file with the ES plots.
#' @return A pdf containing the enrichment plots for all gene sets (with more than 2 genes in a gene set) in a database, saved as the value given plotname.
#' @export

plot.ES=function(list.of.plots="",plotname=""){
  filename <- paste0(plotname,".pdf")
  pdf(filename)
  print(marrangeGrob(list.of.plots, nrow = 2,ncol=2))
  dev.off()
  #title.size<-10
}

