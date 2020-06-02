#' GSEA Enrichment Plots
#'
#' Create a pdf of a list of ggplot objects that is obtained from the GSEA analysis. Outputs 4 graphs per page and saves to your working directory
#' @param list.of.plots A list of ggplot items
#' @param plotname filename to save the plot as
#' @export

plot.ES=function(list.of.plots="",plotname=""){
  filename <- paste0(plotname,".pdf")
  pdf(filename)
  print(marrangeGrob(list.of.plots, nrow = 2,ncol=2))
  dev.off()
  #title.size<-10
}

