% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Run.GSEA.R
\name{GSEAplots}
\alias{GSEAplots}
\title{Execute GSEA Analysis}
\usage{
GSEAplots(
  input.ds.name = "",
  input.cls.name = "",
  gene.set.input = "",
  doc.string = "",
  nperm = 1000,
  bar_percent = 0.1,
  gap_percent = 0.1,
  under_percent = 0.1,
  upper_percent = 0.1,
  color_line = "black",
  color_tick = "black",
  abs.val = F
)
}
\arguments{
\item{input.ds.name}{Name of the gct expression file.}

\item{input.cls.name}{Name of the cls phenotype file.}

\item{gene.set.input}{Name of the geneset file.}

\item{doc.string}{Name of the output folder for analysis and for naming output files.}

\item{nperm}{number of permutations.}

\item{bar_percent}{proportional height of tick mark to window.}

\item{gap_percent}{proportional height between minimum enrichment score and top of tick mark to the window.}

\item{under_percent}{proportional height of white space under tick marks to the window size.}

\item{upper_percent}{proportional height of white space over enrichment graph to the window size.}

\item{color_line}{color of enrichment score line in plot pdf.}

\item{color_tick}{color to tick marks on plot.}

\item{abs.val}{Default is false. Determines whether genes are ranked according to signal to noise or absolute value of signal to noise (when abs.val=T).}
}
\value{
pp - A list pp which includes : plots, gene.set.reference.matrix, gene.set.leading, report1, report2, ES.

plots - A list of ggplot objects, this can act as an input to the function plot.ES() to output a pdf.

gene.set.reference.matrix - A list of each gene set within a gene set database and the gene symbols corresponding to each set.

gene.set.leading - a similar structure to gene.set.reference.matrix but only contains the gene symbols within each gene set that are part of the leading edge set.

report1 - Summary of GSEA analysis data for the first phenotype.

report2 - Summary of GSEA analysis data for the second phenotype.

ES - This object contains the enrichment scores and enrichment tags used to create the plots described earlier. The user can use this information to customize plots as they wish.
}
\description{
Analyzes genetic expression data and determines whether defined gene sets show statistically significant differences with respect to two phenotypes.
}
\details{
GSEA analysis is computed using the Broad Institute's R source code. Genes are ranked according to signal to noise ratio (difference in means/sum of standard deviations for the two phenotypes).
}
\examples{
 pp = GSEAplots(input.ds.name=expr.input, input.cls.name=pheno.input,
      gene.set.input=gene.set.input, doc.string="GSEA_plots", nperm=1000,
      abs.val=F, bar_percent=0.1, gap_percent=0.1, under_percent=0.1,
      upper_percent=0.1, color_line="black", color_tick="black")

 custom_results= GSEAplots(input.ds.name=expr.input, input.cls.name=pheno.input,
      gene.set.input=gene.set.input, doc.string="custom_results", nperm=1000,
      bar_percent=0.1, gap_percent=0.1, under_percent=0.1, upper_percent=0.1,
      color_line="black", color_tick="black", abs.val=F)
}
\references{
Subramanian, Tamayo, et al. (2005), PNAS 102, 15545-15550, http://www.broad.mit.edu/gsea/
}
