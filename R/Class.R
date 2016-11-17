### set class
#' An S4 class that pairs information, probe information and gene information.
#' @slot pairInfo A data.frame
#' @slot probeInfo A GRanges object.
#' @slot geneInfo A GRanges object.
#' @exportClass Pair
setClass("Pair",
         representation = representation(pairInfo="data.frame",
                                         probeInfo="GRanges",geneInfo="GRanges"),
         prototype = list(pairInfo=data.frame(),
                          probeInfo=NULL,geneInfo=NULL),
         validity=function(object){
           message("~~~ Pair: inspector ~~~")
           if(!is.null(probeInfo)){
             if(any(!rownames(object@pairInfo$probe) %in% 
                      as.character(object@probeInfo$name)))
               warning("[Pair: validation] probeInfo doesn't contain information
                       for all probes.")
           }else if(!is.null(geneInfo)){
             if(any(!rownames(object@exp) %in% 
                      as.character(object@geneInfo$GENEID)))
               warning("[Pair: validation] geneInfo doesn't contain 
                       information for all genes.")
           }
           return(TRUE)
         }
)



