# NearGenes
# @param Target A charactor which is name of TRange or one of rownames of TBed.
# @param Gene A GRange object contains coordinates of promoters for human genome.
# @param geneNum A number determine how many gene will be collected from each 
# side of target (number shoule be even).
# @param TRange A GRange object contains coordinate of targets.
# @return A data frame of nearby genes and information: genes' IDs, genes' symbols, 
# distance with target and side to which the gene locate to the target.
#'@importFrom GenomicRanges strand<-
NearGenes <- function (Target=NULL,
                       Gene=NULL,
                       geneNum=20,
                       TRange=NULL){
  
  if(is.null(Gene) | is.null(Target)){
    stop ("Target and Genes should both be defined")
  }
  if(is.null(TRange)){
    stop( "TRange must be defined")
  }else{
    # Just to be sure we have only one probe. To be removed ?
    regionInfo <- TRange[as.character(TRange$name) %in% Target]
  }
  GeneIDs <- c()
  Distances <- c()
  strand(Gene) <- "*"
  # We will get only genes in the same same chromossome
  Gene <- Gene[as.character(seqnames(Gene)) %in% as.character(seqnames(regionInfo))]
  if(length(Gene)==0){
    warning(paste0(Target," don't have any nearby gene in the given gene list."))
    Final <- NA
  }else{
    Gene <- sort(Gene)
    index <- follow(regionInfo,Gene)
    #left side
    Leftlimit <- geneNum/2
    Rightlimit <- geneNum/2
    n <- 1
    if(is.na(index)){
      index<- 0
      Leftlimit <- 0
      Left <- c()
    }else if(index==1){
      Left <- index
      Leftlimit <- length(Left)
    }else{
      Left <- index
      while(length(Left) < Leftlimit){
        if(!as.character(Gene$GENEID[index-n])%in%as.character(Gene$GENEID[Left]))
          Left <- c((index-n),Left) 
        if((index-n)==1){
          Leftlimit <- length(Left)
        } 
        n <- n+1
      }
    }
    
    Right <- c()
    n <- 1
    if(index==length(Gene) || 
       all(unique(Gene$GENEID[(index+1):length(Gene)]) %in% as.character(Gene$GENEID[index]))){
      Rightlimit <- length(Right)
    }else{
      while(length(Right) < Rightlimit){
        if(!as.character(Gene$GENEID[index+n])%in% as.character(Gene$GENEID[c(Right,Left)])) 
          Right <- c(Right,(index+n))
        
        if(index+n==length(Gene)){
          Rightlimit <- length(Right)
        } else{
          n <- n+1
        }     
      }
    }
    
    if(Rightlimit < geneNum/2){
      n <- 1
      if(Left[1]-n > 0){
        while((length(Left)+length(Right)) < geneNum){
          if(!as.character(Gene$GENEID[Left[1]-n])%in%as.character(Gene$GENEID[c(Left,Right)])) 
            Left <- c((Left[1]-n),Left) 
          n <- n+1
        }
      }  
    }
    
    if(Leftlimit < geneNum/2){
      n <- 1
      m <- length(Right)
      if(Right[m]+n < length(Gene)+1)
        while((length(Left)+length(Right)) < geneNum){
          if(!as.character(Gene$GENEID[Right[m]+n])%in%as.character(Gene$GENEID[c(Left,Right)])) 
            Right <- c(Right,(Right[m]+n)) 
          n <- n+1
        } 
    }
    Whole <- c(Left,Right)
    GeneIDs <- Gene$GENEID[Whole]
    Symbols <- Gene$SYMBOL[Whole]
    Distances <-  suppressWarnings(distance(Gene[Whole],regionInfo))
    if(Rightlimit < 1){
      Sides <- paste0("L",length(Left):1)
    }else if( Leftlimit < 1){
      Sides <- paste0("R",1:length(Right))
    }else{
      Sides <- c(paste0("L",length(Left):1),paste0("R",1:length(Right)))
    }
    Final <- data.frame(Target=rep(Target,length(GeneIDs)),GeneID=GeneIDs,
                        Symbol=Symbols,Distance=Distances, Side=Sides, 
                        stringsAsFactors = FALSE)
  }
  return(Final)
}

#' GetNearGenes to collect nearby genes for one locus.
#' @description 
#' GetNearGenes is a function to collect equal number of gene on each side of one locus.
#' @param geneAnnot A GRange object. Contains coordinates of promoters for 
#' human genome.
#' @param geneNum A number determines how many gene will be collected totally. 
#' Then the number devided by 2 is the number of genes collected from 
#' each side of targets (number shoule be even) Default to 20.
#' @param TRange A GRange object. Contains coordinates of a list of targets loci.
#' @param cores A interger which defines the number of cores to be used in parallel process.
#'  Default is 1: no parallel process.
#' @return A data frame of nearby genes and information: genes' IDs, genes' symbols, 
#' distance with target and side to which the gene locate to the target.
#' @export
#' @importFrom GenomicRanges strand follow distance
#' @importFrom pbapply pbsapply
#' @importFrom plyr alply
#' @importFrom doParallel registerDoParallel
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' geneAnnot <- txs(TSS=list(upstream=0, downstream=0))
#' probe <- GRanges(seqnames = c("chr1","chr2"), 
#' range=IRanges(start = c(16058489,236417627), end= c(16058489,236417627)), 
#' name= c("cg18108049","cg17125141"))
#' NearbyGenes <- GetNearGenes(geneNum=20,geneAnnot=geneAnnot,TRange=probe)
GetNearGenes <- function(geneAnnot=NULL,TRange=NULL,geneNum=20,
                         cores=1){
  
  if(is.null(TRange)){
    stop(" TRange must be defined")
  }
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  
  NearGenes <- alply(.data = as.data.frame(TRange), .margins = 1,
                     .fun = function(x) {
                       NearGenes( Target = as.character(x$name),
                                  geneNum = geneNum,
                                  Gene = geneAnnot,
                                  TRange = TRange)},
                     .progress = "text", .parallel = parallel
  )
  names(NearGenes) <- as.data.frame(TRange)$name
  return(NearGenes)
}



