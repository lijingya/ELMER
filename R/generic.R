

#' getMeth
#' @param object MEE.data object
#' @param probe A vector of probes' name. When specified, DNA methylation only 
#' for these probes will be output.
#' @param ID A vector of sample ID. When specified, DNA methylation only for 
#' these samples will be output.
#' @return DNA methylation values.
#' @exportMethod getMeth
#' @docType methods
setGeneric(name="getMeth",
           def=function(object,probe,ID){standardGeneric("getMeth")})



#' getExp
#' @param object MEE.data object
#' @param geneID A vector of genes' id. When specified, gene expression only 
#' for these genes will be output.
#' @param ID A vector of sample ID. When specified, gene expression only for 
#' these samples will be output.
#' @return Gene expression values.
#' @exportMethod getExp
#' @docType methods
setGeneric(name="getExp",
           def=function(object,geneID,ID){standardGeneric("getExp")})


#' getSample
#' @param object MEE.data object
#' @param ID A vector of sample ID. When specified, sample informtion only for 
#' these samples will be output.
#' @param cols A vector of columns names of Sample slots of MEE.data object.
#' @return Sample information.
#' @exportMethod getSample
#' @docType methods
setGeneric(name="getSample",
           def=function(object,ID,cols){standardGeneric("getSample")})

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

#' getProbeInfo
#' @param object MEE.data or Pair object
#' @param probe A vector of probes' name. When specified, only the these probes' 
#' coordinate will be output.
#' @param chr A vector of chromosome such chr1, chr2. When specified, only the 
#' probeInfo locating on these chromosome will be output.
#' @param range A GRanges object. When specified, only the probeInfo locating 
#' within these loci will be output.
#' @return Probe information such as names, coordinates.
#' @exportMethod getProbeInfo
#' @docType methods
setGeneric(name="getProbeInfo",
           def=function(object,chr,probe,range){standardGeneric("getProbeInfo")})



#' getGeneInfo
#' @param object MEE.data or Pair object
#' @param geneID A vector of genes' id. When specified, only the these genes' 
#' coordinate will be output.
#' @param symbol A vector of genes' symbols . When specified, only the these genes' 
#' coordinate will be output.
#' @param range A GRanges object. When specified, only the geneInfo locating within
#'  these loci will be output.
#' @return Gene annotation information such as gene id, symbol and coordinates.
#' @exportMethod getGeneInfo
#' @docType methods
setGeneric(name="getGeneInfo",
           def=function(object,geneID,symbol,range){standardGeneric("getGeneInfo")})
