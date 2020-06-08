% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_phenoinput.R
\name{create_phenoinput}
\alias{create_phenoinput}
\title{Create phenoinput}
\usage{
create_phenoinput(ann = "")
}
\arguments{
\item{ann}{A list of the phenotypes for a sample. There should be two phenotypes.}
}
\value{
A list that is ready to be used in the GSEA.plot function
}
\description{
This function takes two inputs, the phenotype annotations and the order.
}