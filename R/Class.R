### set class
#' MEE.data
#' An S4 class contains 5 slots: methylation, expression, sample information, probe 
#' information and gene information. MEE.data class are the main input for main functions.
#' @slot meth A matrix of DNA methylation. Each row is one probe and each 
#' column is one sample
#' @slot exp A matrix of expression. Each row is one gene and each 
#' column is one sample
#' @slot sample A data.frame contains sample information
#' @slot probeInfo A GRange object contains probe information
#' @slot geneInfo A GRange object contains gene information
#' @import methods GenomicRanges
#' @exportClass MEE.data
setClass("MEE.data",
         representation = representation(meth="matrix",exp="matrix",
                                         sample="data.frame",probeInfo="GRanges"
                                         ,geneInfo="GRanges"),
         validity=function(object){
           message("~~~ MEE.data: inspector ~~~")
           if(!is.null(probeInfo)){
             if(any(!rownames(object@meth) %in% 
                      as.character(object@probeInfo$name)))
               warning("[MEE.data: validation] probeInfo doesn't contain all
                       the information for the probes in the meth array.")
           }else if(!is.null(geneInfo)){
             if(sum(rownames(object@exp) %in% 
                      as.character(object@geneInfo$GENEID))==0)
               stop("[MEE.data: validation] GENEID in geneInfo isn't consistent
                    with expression matrix rownames")
             if(any(!rownames(object@exp) %in% as.character(object@geneInfo$GENEID)))
               warning("[MEE.data: validation] geneInfo doesn't contain all the
                       information for the gene in the expression matrix.")
           }
           return(TRUE)
         }
)


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



