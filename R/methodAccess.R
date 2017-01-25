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


#' getSymbol
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#' @param geneID A character which is the geneID
#' @return The gene symbol 
#' @export
getSymbol <- function(data, geneID){
  gene.info <- values(getExp(data))
  gene <- gene.info[gene.info$ensembl_gene_id %in% geneID,"external_gene_name"]
  return(gene)
}

