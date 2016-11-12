#' getPair to extract pairInfo slot from Pair object.
#' @description 
#' getPair is a function to easily extract pairInfo out of a Pair object.
#' By specifying geneID or probe, geneInfo for the defined genes (geneID ) and probes (probe)
#' will be extracted out of Pair object. 
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

