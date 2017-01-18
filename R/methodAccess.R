#' @title Get DNA methylation object from MAE
#' @description Get DNA methylation object from MAE
#' @importFrom MultiAssayExperiment experiments
#' @export
getMet <- function(data) {
  return(experiments(data)[["DNA methylation"]])
}

#' @title Get DNA methylation object samples from MAE
#' @description Get DNA methylation object samples from MAE
#' @importFrom MultiAssayExperiment sampleMap
#' @export
getMetSamples <- function(data){
  return(sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"])
}

#' @title Get Gene expression object samples from MAE
#' @description Get Gene expression object samples from MAE
#' @importFrom MultiAssayExperiment sampleMap
#' @export
getExpSamples <- function(data){
  return(sampleMap(data)[sampleMap(data)$assay == "Gene expression","primary"])
}

#' @title Get Gene expression object from MAE
#' @description Get Gene expression object from MAE
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


#'getSymbol
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMultiAssayExperiment function}}.
#' @param geneID A character which is the geneID
#' @return The gene symbol 
#' @export
getSymbol <- function(data, geneID){
  gene.info <- values(getExp(data))
  gene <- gene.info[gene.info$ensembl_gene_id %in% geneID,"external_gene_name"]
  return(gene)
}
  

###Pair -------------------------------------------------------------------
# initialize
setMethod(f="initialize",signature="Pair",
          definition=function(.Object,pairInfo,probeInfo,geneInfo){
            message("~~~ Pair: initializator ~~~ ")
            # meth needs to be matrix
            if(!missing(pairInfo)) .Object@pairInfo<-pairInfo
            #probeInfo needs to be GRanges
            if(!missing(probeInfo)) .Object@probeInfo <- probeInfo
            #geneInfo needs to be GRanges
            if(!missing(geneInfo)) .Object@geneInfo <- geneInfo
            return(.Object)
          }
)

#show
setMethod (f="show","Pair",
           function(object){
             cat("*** Class Pair, method show *** \n")
             cat("* pair \n"); print(str(object@pairInfo))
             cat("* probeInfo \n"); print(head(object@probeInfo))
             cat("* geneInfo \n"); print(head(object@geneInfo))
             cat("******* End Print (Pair) ******* \n")
           }
)

setMethod (f="summary","Pair",
           function(object){
             cat("*** Class Pair, method summary *** \n")
             cat("* pair \n"); print(str(object@pairInfo))
             cat("* probeInfo \n"); print(head(object@probeInfo))
             cat("* geneInfo \n"); print(head(object@geneInfo))
             cat("******* End Print (MEE.data) ******* \n")
           }
)


#Accessor Getter

#' @rdname getPair
#' @aliases getPair
#' @examples
#' df <- data.frame(Probe=c("cg19403323","cg12213388","cg26607897"),
#' GeneID =c("ID255928","ID84451","ID55811"),
#' Symbol =c("SYT14","KIAA1804","ADCY10"),
#' Pe=c(0.003322259,0.003322259,0.003322259))
#' pair <- fetch.pair(pair = df)
#' Pairs <- getPair(pair, probe = "cg19403323") # get pair information for a probe
#' Pairs <- getPair(pair, geneID = "ID55811") # get pair information for a gene
setMethod(f="getPair",signature="Pair",
          function(object,probe,geneID){
            if(missing(geneID) & missing(probe)){
              pair <- object@pairInfo
            }else if(missing(probe) & !missing(geneID)){
              pair <- object@pairInfo[object@pairInfo$GeneID %in% geneID,]
            }else if(!missing(probe) & missing(geneID)){
              pair <- object@pairInfo[object@pairInfo$Probe %in% probe,]
            }else if(!missing(probe) & !missing(geneID)){
              pair <- object@pairInfo[object@pairInfo$Probe %in% probe & 
                                        object@pairInfo$GeneID %in% geneID,]
            }
            return(pair)
          })

## Construct data structure-----------------------------------------------------

pair.data <- function(pairInfo=NULL,probeInfo=NULL,geneInfo=NULL){
  args <- list(pairInfo=pairInfo,probeInfo=probeInfo,geneInfo=geneInfo)
  args <- args[!unlist(lapply(args,is.null))]
  args$Class <- "Pair"
  Pair <- do.call(new,args)
  return(Pair)
}


