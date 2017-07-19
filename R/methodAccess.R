#' @title Get DNA methylation object from MAE
#' @description Get DNA methylation object from MAE
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#' @importFrom MultiAssayExperiment experiments
#' @export
getMet <- function(data) {
  return(experiments(data)[["DNA methylation"]])
}

#' @title Get DNA methylation object samples from MAE
#' @description Get DNA methylation object samples from MAE
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#' @importFrom MultiAssayExperiment sampleMap
#' @export
getMetSamples <- function(data){
  return(sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"])
}

#' @title Get Gene expression object samples from MAE
#' @description Get Gene expression object samples from MAE
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#' @importFrom MultiAssayExperiment sampleMap
#' @export
getExpSamples <- function(data){
  return(sampleMap(data)[sampleMap(data)$assay == "Gene expression","primary"])
}

#' @title Get Gene expression object from MAE
#' @description Get Gene expression object from MAE
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#' @importFrom MultiAssayExperiment experiments
#' @export
getExp <- function(data) {
  return(experiments(data)[["Gene expression"]])
}

# Check input
checkData <- function(data){
  if(class(data) == class(MultiAssayExperiment())) {
    if(! "DNA methylation" %in%  names(data) |
       !"Gene expression" %in%  names(data))
      stop("Please the input should be a MultiAssayExperiment with both DNA methylation and Gene expression matrix. See function createMultiAssayExperiment")
  } else {
    stop("Please the input should be a MultiAssayExperiment with both DNA methylation and Gene expression matrix. See function createMultiAssayExperiment")
  }
  # We will need to ensure createMultiAssayExperiment add those fields if they don't exists
  if(!"external_gene_name" %in% colnames(values(data))) stop("Please the input should be a Gene Expression data must have external_gene_name column")
  if(!"ensembl_gene_id" %in% colnames(values(data))) stop("Please the input should be a Gene Expression data must have external_gene_name column")
}

#'getSymbol to report gene symbol from id
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#' @param geneID A character which is the ensembl_gene_id
#' @return The gene symbol for input genes.
#' @export
#' @examples
#' data <- ELMER:::getdata("elmer.data.example")
#' getSymbol(data, geneID="ENSG00000143067")
getSymbol <- function(data,geneID){
  gene <- unique(values(getExp(data))[,c("ensembl_gene_id","external_gene_name")])
  gene <- gene[match(geneID,gene$ensembl_gene_id),"external_gene_name"]
  return(gene)
}

#'getGeneID to report gene id from symbol
#'@importFrom S4Vectors values
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#'@param symbol A vector of characters which are gene symbols 
#'@return The gene ID for these gene symbols
#'@export
#'@examples
#' data <- ELMER:::getdata("elmer.data.example")
#' getGeneID(data, symbol="ZNF697")
getGeneID <- function(data,symbol){
  gene <- unique(values(getExp(data))[,c("ensembl_gene_id","external_gene_name")])
  gene <- gene[match(symbol,gene$external_gene_name),"ensembl_gene_id"]
  return(gene)
}

