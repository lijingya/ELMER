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

#Accessor Getter
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

setGeneric(name="getExp",def=function(object,GeneID,ID){standardGeneric("getExp")})
setMethod(f="getExp",signature="MEE.data",
          function(object,GeneID,ID){
            if(missing(GeneID) & missing(ID)){
              exp <- object@exp
            }else if(missing(GeneID) & !missing(ID)){
              exp <- object@exp[,ID]
            }else if(!missing(GeneID) & missing(ID)){
              exp <- object@exp[GeneID,]
            }else if(!missing(GeneID) & !missing(ID)){
              exp <- object@exp[GeneID,ID]
            }
            return(exp)
          })

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
            if(!missing(pairInfo)){
              pairInfo <- data.frame(probe=as.character(pairInfo$probe),geneID=as.character(pairInfo$geneID),Symbol=as.character(pairInfo$geneName),
                                     Distance = as.numeric(pairInfo$Distance),Side=as.character(pairInfo$Side),Raw.p=as.numeric(pairInfo$Raw.p),Pe=as.numeric(pairInfo$Pe),
                                     stringsAsFactors=F)
              .Object@pairInfo<-pairInfo
            }
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
             cat("******* End Print (MEE.data) ******* \n")
           }
)

#Accessor Getter
setGeneric(name="getPair",def=function(object,pair,cols){standardGeneric("getPair")})
setMethod(f="getPair",signature="Pair",
          function(object,pair,cols){
            if(missing(pair) & missing(cols)){
              pair <- object@pairInfo
            }else if(missing(pair) & !missing(cols)){
              pair <- object@pairInfo[,cols]
            }else if(!missing(pair) & missing(cols)){
              pair <- object@pairInfo[pair,]
            }else if(!missing(pair) & !missing(cols)){
              pair <- object@pairInfo[pair,cols]
            }
            return(pair)
          })

setGeneric(name="getProbeInfo",def=function(object,chr,probe){standardGeneric("getProbeInfo")})
setMethod(f="getProbeInfo",signature=c("MEE.data","Pair"),
          function(object,chr,probe){
            if(missing(chr) & missing(probe)){
              probeInfo <- object@probeInfo
            }else if(missing(chr) & !missing(probe)){
              probeInfo <- object@probeInfo[Object@probeInfo$name %in% probe]
            }else if(!missing(chr) & missing(probe)){
              probeInfo <- object@probeInfo[seqnames(Object@probeInfo) %in% chr]
            }else if(!missing(chr) & !missing(probe)){
              probeInfo <- object@probeInfo[Object@probeInfo$name %in% probe & seqnames(Object@probeInfo) %in% chr]
            }
            return(probeInfo)
          })

setGeneric(name="getGeneInfo",def=function(object,GeneID,GeneName,Simple){standardGeneric("getGeneInfo")})
setMethod(f="getGeneInfo",signature=c("MEE.data","Pair"),
          function(object,GeneID,GeneName,Simple=FALSE){
            if(missing(GeneID) & missing(GeneName)){
              out <- object@geneInfo
            }else if(missing(GeneID) & !missing(GeneName)){
              if(sample){
                out <- unique(as.character(object@geneInfo$GENEID[object@geneInfo$SYMBOL %in% GeneName]))
              }else{
                out <- object@geneInfo[object@geneInfo$SYMBOL %in% GeneName]
              }
            }else if(!missing(GeneID) & missing(GeneName)){
              if(Simple){
                out <- unique(as.character(object@geneInfo$SYMBOL[object@geneInfo$GENEID %in% GeneID]))
              }else{
                out <- object@geneInfo[object@geneInfo$GENEID %in% GeneID]
              }
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


