#' Collect nearby gene for one locus.
#' @param Target A charactor which is name of TRange or one of rownames of TBed.
#' @param Gene A GRange object contains coordinates of promoters for human genome.
#' @param SampleSize A number determine how many gene will be collected from each side of target (number shoule be even).
#' @param Tbed A bed format data.frame object contains coordinate of targets.
#' @param TRange A GRange object contains coordinate of targets.
#' @return A data frame of nearby genes and information: genes' IDs, genes' symbols, distance with target and side to which the gene locate to the target.
 
.NearGenes <- function (Target=NULL,Gene=NULL,SampleSize=20,TBed=NULL,TRange=NULL){
  if(is.null(Gene) | is.null(Target)){
    stop ("Target and Genes should both be defined")
  }
  source("/export/uec-gs1/laird/users/lijingya/software/scripts/R/ReadFile.R")
  require("GenomicRanges")
  message(Target)
  # form the Gene GRange
  # Gene <- GRanges(seqnames = TBed[Target,1], ranges = IRanges (start =as.numeric(TBed[Target,2]), end = as.numeric(TBed[Target,3])), strand=TBed[Target,6])
  if(is.null(TRange)){
    if(is.null(TBed)){
      stop( "Either TBed or TRange must be defined")
    }
    regionInfo <- GRanges(seqnames = TBed[Target,1], ranges = IRanges (start =as.numeric(TBed[Target,2]), end = as.numeric(TBed[Target,3])))
  }else{
    regionInfo <- TRange[TRange$name %in% Target,]
  }
  GeneIDs <-c()
  Distances <- c()
  strand(Gene) <- "*"
  Gene <- Gene[as.character(seqnames(Gene)) %in% as.character(seqnames(regionInfo))]
  if(length(Gene)==0){
    warning(paste0(Target," don't have any nearby gene in the given gene list."))
    Final <- NA
  }else{
    Gene <- sort(Gene)
    #  Gene_subset <- Gene[as.character(seqnames(Gene)) %in% as.character(seqnames(regionInfo))]
    index <- follow(regionInfo,Gene)
    #left side
    Leftlimit <- SampleSize/2
    Rightlimit <- SampleSize/2
    Left <- index
    n <- 1
    if(index==1){
      Leftlimit <- length(Left)
    }else{
      while(length(Left) < Leftlimit){
        if(!as.character(Gene$GENEID[index-n])%in%as.character(Gene$GENEID[Left])) Left <- c((index-n),Left) 
        if((index-n)==1){
          Leftlimit <- length(Left)
        } 
        n <- n+1
      }
    }
    
    Right <- c()
    n <- 1
    if(index==length(Gene) || all(unique(Gene$GENEID[(index+1):length(Gene)]) %in% as.character(Gene$GENEID[index]))){
      Rightlimit <- length(Right)
    }else{
      while(length(Right) < Rightlimit){
        if(!as.character(Gene$GENEID[index+n])%in% as.character(Gene$GENEID[c(Right,Left)])) Right <- c(Right,(index+n))
        
        if(index+n==length(Gene)){
          Rightlimit <- length(Right)
        } else{
          n <- n+1
        }     
      }
    }
    
    if(Rightlimit < SampleSize/2){
      n <- 1
      if(Left[1]-n > 0){
        while((length(Left)+length(Right)) < SampleSize){
          if(!as.character(Gene$GENEID[Left[1]-n])%in%as.character(Gene$GENEID[c(Left,Right)])) Left <- c((Left[1]-n),Left) 
          n <- n+1
        }
      }  
    }
    
    if(Leftlimit < SampleSize/2){
      n <- 1
      m <- length(Right)
      if(Right[m]+n < length(Gene)+1)
        while((length(Left)+length(Right)) < SampleSize){
          if(!as.character(Gene$GENEID[Right[m]+n])%in%as.character(Gene$GENEID[c(Left,Right)])) Right <- c(Right,(Right[m]+n)) 
          n <- n+1
        } 
    }
    print(Left)
    print(Right)
    Whole <- c(Left,Right)
    GeneIDs <- Gene$GENEID[Whole]
    Symbols <- Gene$SYMBOL[Whole]
    Distances <-  distance(Gene[Whole],regionInfo)
    if(Rightlimit < 1){
      Sides <- paste0("L",length(Left):1)
    }else if( Leftlimit < 1){
      Sides <- paste0("R",1:length(Right))
    }else{
      Sides <- c(paste0("L",length(Left):1),paste0("R",1:length(Right)))
    }
    Final <- cbind(Target=rep(Target,length(GeneIDs)),GeneID=GeneIDs,Symbol=Symbols,Distance=Distances, Sides=Sides)
    names(Final) <- NULL
  }
  return(Final)
}

#' Collect nearby gene for one locus.
#' @param Gene A GRange object contains coordinates of promoters for human genome.
#' @param SampleSize A number determine how many gene will be collected from each side of target (number shoule be even) Default to 20.
#' @param Tbed A bed format data.frame object contains coordinate of a list targets.
#' @param TRange A GRange object contains coordinate of a list targets.
#' @param cores A number to specific how many cores to use to compute. Default to detectCores()/2.
#' @return A data frame of nearby genes and information: genes' IDs, genes' symbols, distance with target and side to which the gene locate to the target.

GetNearGenes <- function(SampleSize=20,Gene=NULL,TBed=NULL,TRange=NULL,cores=NULL){
  require("snow")
  if(is.null(cores)) cores <- detectCores()/2   
  cl <- makeCluster(cores,type="SOCK")
  #  options('mc.cores'=detectCores()/2)
  if(is.null(TRange)){
    if(is.null(TBed)){
      stop( "Either TBed or TRange must be defined")
    }
    rownames(TBed) <- TBed[,4]
    out <- parSapplyLB(cl,rownames(TBed),.NearGenes,SampleSize=SampleSize,Gene=Gene,TBed=TBed,simplify=FALSE)
  }else{
    out <- parSapplyLB(cl,as.character(TRange$name),.NearGenes,SampleSize=SampleSize,Gene=Gene,TRange=TRange,simplify=FALSE)
  }
  stopCluster(cl)
  return(out)
}

