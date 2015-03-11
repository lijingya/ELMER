###MEE.data ----------------------------------------------------
# # initialize 
setMethod(f="initialize",signature="MEE.data",
          definition=function(.Object,meth,exp,sample,probeInfo,geneInfo){
            cat("~~~ MEE.data: initializator ~~~ \n")
            # meth needs to be matrix
            if(!missing(meth)) .Object@meth<-meth
            #exp needs to be matrix
            if(!missing(exp)) .Object@exp<-exp
            #sample needs to be data.frame
            if(!missing(sample)) .Object@sample <- sample
            #probeInfo needs to be GRanges
            if(!missing(probeInfo)) .Object@probeInfo <- probeInfo
            #geneInfo needs to be GRanges
            if(!missing(geneInfo)) .Object@geneInfo <- geneInfo
            return(.Object)
          }
)

# show
setMethod (f="show","MEE.data",
           function(object){
             cat("*** Class MEE.data, method show *** \n")
             cat("* meth \n"); print(str(object@meth))
             cat("* exp \n"); print(str(object@exp))
             cat("* sample \n"); print(str(object@sample))
             cat("* probeInfo \n"); print(object@probeInfo)
             cat("* geneInfo \n"); print(object@geneInfo)
             cat("******* End Print (MEE.data) ******* \n")
           }
)

setMethod (f="summary","MEE.data",
           function(object){
             cat("*** Class MEE.data, method summary *** \n")
             cat("* meth \n"); print(str(object@meth))
             cat("* exp \n"); print(str(object@exp))
             cat("* sample \n"); print(str(object@sample))
             cat("* probeInfo \n"); print(object@probeInfo)
             cat("* geneInfo \n"); print(object@geneInfo)
             cat("******* End Print (MEE.data) ******* \n")
           }
)

#Accessor Getter
#' getMeth
#' @param object MEE.data object
#' @param probe A vector of probes' name. When specified, DNA methylation only for these probes will be output.
#' @param ID A vector of sample ID. When specified, DNA methylation only for these samples will be output.
#' @exportMethod getMeth
setGeneric(name="getMeth",def=function(object,probe,ID){standardGeneric("getMeth")})
setMethod(f="getMeth",signature="MEE.data",
          definition=function(object,probe,ID){
            if(missing(probe) & missing(ID)){
              meth <- object@meth
            }else if(missing(probe) & !missing(ID)){
              meth <- object@meth[,ID]
            }else if(!missing(probe) & missing(ID)){
              meth <- object@meth[probe,]
            }else if(!missing(probe) & !missing(ID)){
              meth <- object@meth[probe,ID]
            }
            return(meth)
          })

#' getExp
#' @param object MEE.data object
#' @param geneID A vector of genes' id. When specified, gene expression only for these genes will be output.
#' @param ID A vector of sample ID. When specified, gene expression only for these samples will be output.
#' @exportMethod getExp

setGeneric(name="getExp",def=function(object,geneID,ID){standardGeneric("getExp")})
setMethod(f="getExp",signature="MEE.data",
          function(object,geneID,ID){
            if(missing(geneID) & missing(ID)){
              exp <- object@exp
            }else if(missing(geneID) & !missing(ID)){
              exp <- object@exp[,ID]
            }else if(!missing(geneID) & missing(ID)){
              exp <- object@exp[geneID,]
            }else if(!missing(geneID) & !missing(ID)){
              exp <- object@exp[geneID,ID]
            }
            return(exp)
          })

#' getSample
#' @param object MEE.data object
#' @param ID A vector of sample ID. When specified, sample informtion only for these samples will be output.
#' @param cols A vector of columns names of Sample slots of MEE.data object.
#' @exportMethod getSample
#' 
setGeneric(name="getSample",def=function(object,ID,cols){standardGeneric("getSample")})
setMethod(f="getSample",signature="MEE.data",
          function(object,ID,cols){
            if(missing(ID) & missing(cols)){
              sample <- object@sample
            }else if(missing(ID) & !missing(cols)){
              sample <- object@sample[,cols]
            }else if(!missing(ID) & missing(cols)){
              sample <- object@sample[ID,]
            }else if(!missing(ID) & !missing(cols)){
              sample <- object@sample[ID,cols]
            }
            return(sample)
          })





###Pair -------------------------------------------------------------------
# initialize
setMethod(f="initialize",signature="Pair",
          definition=function(.Object,pairInfo,probeInfo,geneInfo){
            cat("~~~ Pair: initializator ~~~ \n")
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
#' getPair
#' @param object Pair object
#' @param geneID A vector of genes' id. When specified, only the pair containing these genes will be output.
#' @param probe A vector of probes' name. When specified, only the pair containing these probes will be output.
#' @exportMethod getPair

setGeneric(name="getPair",def=function(object,probe,geneID){standardGeneric("getPair")})
setMethod(f="getPair",signature="Pair",
          function(object,probe,geneID){
            if(missing(geneID) & missing(probe)){
              pair <- object@pairInfo
            }else if(missing(probe) & !missing(geneID)){
              pair <- object@pairInfo[object@pairInfo$GeneID %in% geneID,]
            }else if(!missing(probe) & missing(geneID)){
              pair <- object@pairInfo[object@pairInfo$Probe %in% probe,]
            }else if(!missing(probe) & !missing(geneID)){
              pair <- object@pairInfo[object@pairInfo$Probe %in% probe & object@pairInfo$GeneID %in% geneID,]
            }
            return(pair)
          })

#' getProbeInfo
#' @param object MEE.data or Pair object
#' @param probe A vector of probes' name. When specified, only the these probes' coordinate will be output.
#' @param chr A vector of chromosome such chr1, chr2. When specified, only the probeInfo locating on these chromosome will be output.
#' @param range A GRanges object. When specified, only the probeInfo locating within these loci will be output.
#' @exportMethod getProbeInfo
setGeneric(name="getProbeInfo",def=function(object,chr,probe,range){standardGeneric("getProbeInfo")})
setMethod(f="getProbeInfo",signature="ANY",
          function(object,chr,probe,range){
            if(missing(chr) & missing(probe)){
              probeInfo <- object@probeInfo
            }else if(missing(chr) & !missing(probe)){
              probeInfo <- object@probeInfo[object@probeInfo$name %in% probe]
            }else if(!missing(chr) & missing(probe)){
              probeInfo <- object@probeInfo[seqnames(object@probeInfo) %in% chr]
            }else if(!missing(chr) & !missing(probe)){
              probeInfo <- object@probeInfo[object@probeInfo$name %in% probe & seqnames(object@probeInfo) %in% chr]
            }
            if(!missing(range)){
              over <- findOverlaps(probeInfo,range)
              probeInfo <- probeInfo[unique(queryHits(over))]
            }
            return(probeInfo)
          })

#' getProbeInfo
#' @param object MEE.data or Pair object
#' @param geneID A vector of genes' id. When specified, only the these genes' coordinate will be output.
#' @param symbol A vector of genes' symbols . When specified, only the these genes' coordinate will be output.
#' @param range A GRanges object. When specified, only the geneInfo locating within these loci will be output.
#' @exportMethod getGeneInfo
setGeneric(name="getGeneInfo",def=function(object,geneID,symbol,range){standardGeneric("getGeneInfo")})
setMethod(f="getGeneInfo",signature="ANY",
          function(object,geneID,symbol,range){
            if(missing(geneID) & missing(symbol)){
              out <- object@geneInfo
            }else if(missing(geneID) & !missing(symbol)){
              out <- object@geneInfo[object@geneInfo$SYMBOL %in% symbol]
            }else if(!missing(geneID) & missing(symbol)){
              out <- object@geneInfo[object@geneInfo$GENEID %in% geneID]
            }
            if(!missing(range)){
              over <- findOverlaps(out,range)
              out <- out[unique(queryHits(over))]
            }
            return(out)
          })

## Construct data structure-----------------------------------------------------

mee.data <- function(meth=NULL,exp=NULL,sample=NULL,probeInfo=NULL,geneInfo=NULL){
  args <- list(meth=meth,exp=exp,sample=sample,probeInfo=probeInfo,geneInfo=geneInfo)
  args <- args[!unlist(lapply(args,is.null))]
  args$Class <- "MEE.data"
  MethExp <- do.call(new,args)
  return(MethExp)
}

pair.data <- function(pairInfo=NULL,probeInfo=NULL,geneInfo=NULL){
  args <- list(pairInfo=pairInfo,probeInfo=probeInfo,geneInfo=geneInfo)
  args <- args[!unlist(lapply(args,is.null))]
  args$Class <- "Pair"
  Pair <- do.call(new,args)
  return(Pair)
}


