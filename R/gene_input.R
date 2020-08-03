#' View Available Geneset Databases
#'
#' The function database_key returns the name of the rda file the user may be interested in given the common name of the geneset database.
#' The function may also return the name of all databases if the input is "all".
#'
#' @param database_name The input for this paramter needs to be "all", or the description of the database of interest.
#' 
#' @return If the input for database_name is "all", the function will return the descriptons for all databases. If the user inputs the descripton of a database, the function will return the file name of the database of interest.
#' @examples
#' database_key("all")
#' 
#' database_key("hallmark all")
#' 
#' database_key("motif transcription factor targets")
#' @export


database_key=function(database_name=""){
  data(key)
  if (database_name=="all"){
    return(as.character(key[,2]))
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
#' This function returns the names of the sets within a database given the file name for the database.
#'
#' @param database_file Name of the rda file for the geneset database within the package.
#' @return The names of the genesets within the database.
#' @examples
#' get_genesets(hallmark.gs)
#' 
#' get_genesets(C3.tft.gs)
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


#' Format New Geneset Database
#'
#' The create_geneset_db function allows users to create a new database. Given a list wich meets certain specifications, this function will reformat the data in the list to be saved into a database form compatible with GSEAplot functions.
#' @param database A list, with the name of the list being the names of the genesets and elements of the list as the description/source and gene symbols.
#' @return Database, saved as formatted_db, which may be saved to the work environment and an used in later GSEA analyses.
#' @examples 
#' create_geneset_db(custom_db)
#' 
#' create_geneset_db(modified_db)
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
#' The get_genesymbols function returns the gene symbols within each gene set given the file name of the database of interest.
#'
#' @param database_file Name of the rda file for the geneset database within the package.
#' @return The gene symbols within each gene set in a database.
#' @examples 
#' get_genesymbols(hallmark.gs)
#' 
#' get_genesymbols(C3.tft.gs)
#' 
#' get_genesymbols(custom_db)
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

#' Add to Existing Gene Set Database
#'
#' The function add_to_database returns a database comprised of an existing database and new data, which are paramters provded by the user. The dataset that results from this function will be formatted to be compatible with GSEA functions.
#' 
#'
#' @param database Formatted existing database.
#' @param addition Gene sets to be added to existing database.
#' @return A formatted database with new sets included.
#' @examples 
#' add_to_database(database=hallmark.gs,addition=new_geneset)
#' 
#' add_to_database(database=custom_db,addition=new_geneset)
#' @export

add_to_database=function(database="",addition=""){
  original_length=length(database)
  for (ii in 1: nrow(addition)){
    temp=paste(addition[ii,], collapse="\t")
    database[[original_length+ii]]=temp
  }
  return(database)
}

