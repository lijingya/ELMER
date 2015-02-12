### set class
#' MEE.data
#' An S4 class that methylation, expression, sample information, probe information and gene information.
#' @slot meth A matrix of DNA methylation. Each row is one probe and each column is one sample
#' @slot exp A matrix of expression. Each row is one gene and each column is one sample
#' @slot sample A data.frame contains sample information
#' @slot probeInfo A GRange object contains probe information
#' @slot geneInfo A GRange object contains gene information
#' @import GenomicRanges
#' @export
setClass("MEE.data",
         representation = representation(meth="matrix",exp="matrix",sample="data.frame",probeInfo="GRanges",geneInfo="GRanges"),
         validity=function(object){
           cat("~~~ MEE.data: inspector ~~~\n")
           if(!is.null(probeInfo)){
             if(any(!rownames(object@meth) %in% as.character(object@probeInfo$name)))
               warning("[MEE.data: validation] probeInfo doesn't contain all the information for the probes in the meth array.")
           }else if(!is.null(geneInfo)){
             if(any(!rownames(object@exp) %in% as.character(object@geneInfo$GENEID)))
               warning("[MEE.data: validation] geneInfo doesn't contain all the information for the gene in the exp array.")
           }
           return(TRUE)
         }
)


### set class
#' An S4 class that pairs information, probe information and gene information.
#' @slot pairInfo A data.frame
#' @slot probeInfo A GRanges object.
#' @slot geneInfo A GRanges object.
#' @export
setClass("Pair",
         representation = representation(pairInfo="data.frame",probeInfo="GRanges",geneInfo="GRanges"),
         prototype = list(pairInfo=data.frame(probe=character(),geneID=character(),Symbol=character(),Distance=numeric(),Side=character(),Raw.p=numeric(),Pe=numeric(),stringsAsFactors=F),
                          probeInfo=NULL,geneInfo=NULL),
         validity=function(object){
           cat("~~~ Pair: inspector ~~~\n")
           if(!is.null(probeInfo)){
             if(any(!rownames(object@pairInfo$probe) %in% as.character(object@probeInfo$name)))
               warning("[Pair: validation] probeInfo doesn't contain information for all probes.")
           }else if(!is.null(geneInfo)){
             if(any(!rownames(object@exp) %in% as.character(object@geneInfo$GENEID)))
               warning("[Pair: validation] geneInfo doesn't contain information for all genes.")
           }
           return(TRUE)
         }
)







###