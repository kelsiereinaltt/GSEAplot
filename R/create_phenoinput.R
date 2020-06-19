#' Create phenoinput
#'
#' The create_phenoinput function translates phenotypes into a format computable by the GSEA algorithm given a list of phenotypes. 
#' 
#'
#' @param ann A list of the phenotypes for a sample. There should be two phenotypes (e.g., female and male).
#' @return A list of the phenotypes converted into a GSEAplot computable format 
#' @export

create_phenoinput=function(ann=""){
  phen=unique(ann)
  first=gsub(phen[1],0,ann)
  second=gsub(phen[2],1,first)
  class.v=as.double(second)

  pheno.input=list(phen=phen, class.v=class.v)
  return(pheno.input)
  }

