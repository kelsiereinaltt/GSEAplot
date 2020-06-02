#' Create phenoinput
#'
#' This function takes two inputs, the phenotype annotations and the order.
#'
#' @param ann A list of the phenotypes for a sample. There should be two phenotypes.
#' @return A list that is ready to be used in the GSEA.plot function
#' @export

create_phenoinput=function(ann=""){
  phen=unique(ann)
  first=gsub(phen[1],0,ann)
  second=gsub(phen[2],1,first)
  class.v=as.double(second)

  pheno.input=list(phen=phen, class.v=class.v)
  return(pheno.input)
  }

