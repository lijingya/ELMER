#' getPair
#' @param object Pair object
#' @param geneID A vector of genes' id. When specified, only the pair containing 
#' these genes will be output.
#' @param probe A vector of probes' name. When specified, only the pair containing 
#' these probes will be output.
#' @return Pair information such as empirical P values, probe and gene ID.
#' @exportMethod getPair
#' @docType methods
setGeneric(name="getPair",
           def=function(object,probe,geneID){standardGeneric("getPair")})

