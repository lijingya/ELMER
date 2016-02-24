###MEE.data ----------------------------------------------------
# # initialize 
setMethod(f="initialize",signature="MEE.data",
          definition=function(.Object,meth,exp,sample,probeInfo,geneInfo){
            message("~~~ MEE.data: initializator ~~~ ")
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
#' @rdname getMeth
#' @aliases getMeth
#' @examples 
#' meth <- matrix(data=c(1:20),ncol=5,dimnames=list(paste0("probe",1:4),paste0("sample",1:5)))
#' mee <- fetch.mee(meth=meth)
#' Meth <- getMeth(mee,probe = "probe1")
#' Meth <- getMeth(mee, ID = c("sample1","sample2"))
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


#' @rdname getExp
#' @aliases getExp
#' @examples
#' exp <- matrix(data=c(101:110),ncol=5,dimnames=list(c("gene1","gene2"),paste0("sample",1:5)))
#' mee <- fetch.mee(exp=exp)
#' Exp <- getExp(mee, geneID = "gene1") ## get gene expression for certain genes
#' Exp <- getExp(mee, ID = c("sample1","sample5")) ## get gene expression for certain samples
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


#' @rdname getSample
#' @aliases getSample
#' @examples
#' SampleInfo <- data.frame(ID=paste0("sample",1:5), 
#' TN=c("Experiment","Experiment","Control","Control","Experiment"))
#' mee <- fetch.mee(sample = SampleInfo)
#' Samples <- getSample(mee,ID = "sample2") ## get sample2's information
#' Samples <- getSample(mee, cols = "TN")  ## get 'TN' information for each samples
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

#' @rdname getProbeInfo
#' @aliases getProbeInfo
#' @import IRanges
#' @examples
#' probeInfo <- GRanges(seqnames = c("chr1","chr1","chr3"), 
#' ranges = IRanges(start = c(1,6,20),end = c(2,7,21)),
#' name=c("cg1","cg2","cg3"))
#' mee <- fetch.mee(probeInfo=probeInfo)
#' Probes <- getProbeInfo(mee,chr="chr1") # get probes which locate on the chr1
#' Probes <- getProbeInfo(mee, probe = "cg1") # get certain probes information
#' Probes <- getProbeInfo(mee, range = GRanges(seqnames="chr1", ranges=IRanges(5,20)))
setMethod(f="getProbeInfo",signature="ANY",
          function(object,chr,probe,range){
            if(missing(chr) & missing(probe)){
              probeInfo <- object@probeInfo
            }else if(missing(chr) & !missing(probe)){
              probeInfo <- object@probeInfo[object@probeInfo$name %in% probe]
            }else if(!missing(chr) & missing(probe)){
              probeInfo <- object@probeInfo[as.vector(seqnames(object@probeInfo)) %in% chr]
            }else if(!missing(chr) & !missing(probe)){
              probeInfo <- object@probeInfo[object@probeInfo$name %in% probe & 
                                              as.vector(seqnames(object@probeInfo)) %in% chr]
            }
            if(!missing(range)){
              over <- findOverlaps(probeInfo,range)
              probeInfo <- probeInfo[unique(queryHits(over))]
            }
            return(probeInfo)
          })

#' @rdname getGeneInfo
#' @aliases getGeneInfo
#'@import methods S4Vectors IRanges GenomicRanges
#'@examples
#'geneInfo <- txs()
#'mee <- fetch.mee(geneInfo=geneInfo)
#'Genes <- getGeneInfo(mee, geneID = "55811")
#'Genes <- getGeneInfo(mee, symbol ="ADCY10")
#'Genes <- getGeneInfo(mee, range = GRanges(seqnames="chr1", ranges=IRanges(1000000,1600000)))
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


