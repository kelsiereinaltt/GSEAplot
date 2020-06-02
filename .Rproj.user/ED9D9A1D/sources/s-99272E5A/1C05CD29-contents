#' View Available Geneset Databases
#'
#' This function takes one input, the common name of the geneset database you are interested in and returns the name of the rda file that can be used in the GSEA analysis function.
#' If the input is "all" then all database common names are returned.
#'
#' @param database_name Either enter "all" or a string containing the name of one of the existing geneset databases
#' @return Either the names of all existing gene set databases or the file name for the specified database
#' @export


database_key=function(database_name=""){
  data(key)
  if (database_name=="all"){
    return(key[,2])
  }
  else {
    index=which(key[,2]==database_name)
    database_of_interest=key[index,1]
    database_of_interest=as.character(database_of_interest)
    return(database_of_interest)
  }
}

#' View Gene Sets within Database
#'
#' This function takes one input, the file name for the database, and returns the names of the sets within the database.
#'
#' @param database_file Name of the rda file for the geneset database within the package
#' @return The names of the genesets within the database
#' @export

get_genesets=function(database_file=""){
temp=database_file
max.Ng <- length(temp)
names <- vector(length = max.Ng, mode = "character")
gs.count <- 1


for (i in 1:max.Ng) {
  gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
  gene.set.name <- gs.line[1]
  names[gs.count] <- gene.set.name
  gs.count <- gs.count + 1
}
return(names)
}

#' Format New Geneset Daatbase
#'
#' When given a list with keys as a set name and the data corresponding to a key, the gene symbols corresponding to that key. This function will reformat the data to be saved into the package and to be a working input to the GSEA.plots function.
#' @param database a list names are the sets within the database and contained within those are genesymbols
#' @return formatted_db a database ready to be saved to the package and used in later GSEA analyses
#' @export
create_geneset_db=function(database=""){
  formatted_db=list()
  names=names(database)

  for (i in 1: length(database)){
    names(database[[i]])=NULL
    q=paste(database[[i]],collapse="\t")
    formatted_db[[i]]=paste(names[[i]],q,sep="\t")
  }

  formatted_db=unlist(formatted_db)
  return(formatted_db)
}

#' Get Gene Symbols
#'
#' This function takes one input, the file name for the database and and returns the gene symbols within that set.
#'
#' @param database_file Name of the rda file for the geneset database within the package
#' @return The gene symbols within that gene set
#' @export

get_genesymbols=function(database_file=""){
  temp=database_file
  max.Ng <- length(temp)
  names <- vector(length = max.Ng, mode = "character")
  db <- list(length=max.Ng,mode="character")
  gs.count <- 1


  for (i in 1:max.Ng) {
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1]
    names[gs.count] <- gene.set.name
    gene.symbols <- gs.line[-1]
    gene.symbols <- gene.symbols[-1]
    db[[gs.count]] <- gene.symbols
    gs.count <- gs.count + 1
  }
  names(db)=names
  return(db)
}

#' Add to Existing Gene set Database
#'
#' This function takes two inputs: the formatted existing database and the gene sets to add.It will return the combined database ready to be used in the GSEA function.
#'
#' @param database the formatted existing database
#' @param addition the gene sets to be added
#' @return a formatted database with new sets included
#' @export

add_to_database=function(database="",addition=""){
  original_length=length(database)
  for (ii in 1: nrow(addition)){
    temp=paste(addition[ii,], collapse="\t")
    database[[original_length+ii]]=temp
  }
  return(database)
}

