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
  # Algorithm:
  # 1) get the follow gene (overlapping genes are diconsidered) to be the first in L1 (index variable)
  #                 probe
  #     -------       O      ------
  #    |      |      ||     |     |
  #    -------       ||     ------
  #  follow gene          precede gene
  # 2) Sort genes (by start) 
  # 3.1) Get 9 genes before index (L1)
  #      If we only have l genes (l < 9) due to end of genomic region, get the l ones and get more 10-l to the right
  # 3.2) Get 10 after index (L1)
  #      If we only have r genes (r < 10) due to end of genomic region, get the r ones and get more 10-r to the left 
  # Where 10 is genum/2
  Gene$GENEID <- Gene$ensembl_gene_id
  if("external_gene_name" %in% colnames(S4Vectors::mcols(Gene))){
    Gene$SYMBOL <- Gene$external_gene_name
  } else  if("external_gene_id" %in% colnames(S4Vectors::mcols(Gene))){
    Gene$SYMBOL <- Gene$external_gene_id
  } else {
    stop("No gene symbol column found (expected external_gene_id or external_gene_name")
  }
  if(is.null(Gene) | is.null(Target)){
    stop ("Target and Genes should both be defined")
  }
  if(is.null(TRange)){
    stop( "TRange must be defined")
  }else{
    # Just to be sure we have only one probe. To be removed ?
    regionInfo <- TRange[names(TRange) %in% Target]
  }
  GeneIDs <- c()
  Distances <- c()
  strand(Gene) <- "*"
  # We will get only genes in the same same chromossome
  Gene <- Gene[as.character(seqnames(Gene)) %in% as.character(seqnames(regionInfo))]
  #print(Gene)
  if(length(Gene)==0){
    warning(paste0(Target," don't have any nearby gene in the given gene list."))
    Final <- NA
  } else {
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
        # If the gene is not in the list already add it, otherwise go to the next
        if(!as.character(Gene$GENEID[index-n]) %in% as.character(Gene$GENEID[Left])) Left <- c((index-n),Left) 
        
        # Is it the first gene? If so there is nothing in the left anymore
        if((index-n)==1) Leftlimit <- length(Left)
        n <- n + 1
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
    } else if( Leftlimit < 1){
      Sides <- paste0("R",1:length(Right))
    } else{
      Sides <- c(paste0("L",length(Left):1),paste0("R",1:length(Right)))
    }
    
    Final <- data.frame(Target=rep(Target,length(GeneIDs)),GeneID=GeneIDs,
                        Symbol=Symbols,Distance=Distances, Side=Sides, 
                        stringsAsFactors = FALSE)
    Final <- Final[order(Final$Side,Final$Distance),]
  }
  return(Final)
}

#' GetNearGenes to collect nearby genes for one locus.
#' @description 
#' GetNearGenes is a function to collect equal number of gene on each side of one locus.
#' It can receite either multi Assay Experiment with both DNA methylation and gene Expression matrix
#' and the names of probes to select nearby genes, or it can receive two granges objects TRange and geneAnnot.
#' @param data A multi Assay Experiment with both DNA methylation and gene Expression objects
#' @param probes Name of probes to get nearby genes (it should be rownames of the DNA methylation 
#' object in the data argument object)
#' @param geneAnnot A GRange object  or Summarized Experiment object that contains coordinates of promoters for 
#' human genome.
#' @param numFlankingGenes A number determines how many gene will be collected totally. 
#' Then the number devided by 2 is the number of genes collected from 
#' each side of targets (number shoule be even) Default to 20.
#' @param TRange A GRange object or Summarized Experiment object that contains coordinates of a list of targets loci.
#' @param cores A interger which defines the number of cores to be used in parallel process.
#'  Default is 1: no parallel process.
#' @return A data frame of nearby genes and information: genes' IDs, genes' symbols, 
#' distance with target and side to which the gene locate to the target.
#' @export
#' @importFrom GenomicRanges strand follow distance 
#' @importFrom plyr alply
#' @importFrom doParallel registerDoParallel
#' @importFrom SummarizedExperiment rowRanges
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' geneAnnot <- getTSS(TSS=list(upstream=0, downstream=0))
#' probe <- GenomicRanges::GRanges(seqnames = c("chr1","chr2"), 
#' range=IRanges::IRanges(start = c(16058489,236417627), end= c(16058489,236417627)), 
#' name= c("cg18108049","cg17125141"))
#' names(probe) <- c("cg18108049","cg17125141")
#' NearbyGenes <- GetNearGenes(numFlankingGenes = 20,geneAnnot=geneAnnot,TRange=probe)
GetNearGenes <- function(data = NULL,
                         probes = NULL,
                         geneAnnot = NULL,
                         TRange = NULL,
                         numFlankingGenes = 20,
                         cores = 1){
  message("Searching for the ", numFlankingGenes, " near genes")
  if(!is.null(data)){
    if(is.null(probes)) stop("Please set the probes argument (names of probes to select nearby genes)")
    TRange <- subset(getMet(data), rownames(getMet(data)) %in% probes)
    geneAnnot <- getExp(data)
  }    
  if(is.null(TRange)){
    stop("TRange must be defined")
  }
  if(is.null(geneAnnot)){
    stop("geneAnnot must be defined")
  }
  
  if(class(TRange) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
    TRange <- rowRanges(TRange)
  }
  if(class(geneAnnot) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
    geneAnnot <- rowRanges(geneAnnot)
  }
  
  parallel <- FALSE
  if (cores > 1){
    if (cores > parallel::detectCores()) cores <- parallel::detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  
  if(is.null(names(TRange))) {
    if(is.null(TRange$name)) stop("No probe names found in TRange")
    names(TRange) <- TRange$name
  }
  
  NearGenes <- alply(.data = as.data.frame(TRange), .margins = 1,
                     .fun = function(x) {
                       NearGenes( Target = rownames(x),
                                  geneNum = numFlankingGenes,
                                  Gene = geneAnnot,
                                  TRange = TRange)},
                     .progress = "text", .parallel = parallel
  )
  names(NearGenes) <- names(TRange)
  if("split_labels" %in% names(attributes(NearGenes))) attributes(NearGenes)["split_labels"] <- NULL 
  if("split_type" %in% names(attributes(NearGenes))) attributes(NearGenes)["split_type"] <- NULL 
  
  if(!is.null(data)) {
    NearGenes$Distance <- NULL
    # For a given probe and gene find nearest TSS
    tss <- getTSS(metadata(data)$genome)
    # To calculate the distance we will get the transcript list
    # Consider only the TSS  (start of the transcript single base)
    tss <- promoters(tss, upstream = 0, downstream = 0)
    message("Update the distance to gene to distance to the nearest TSS of the gene")
    for(probe in names(NearGenes)){
      aux <- NearGenes[[probe]]
      distance <-   plyr::ddply(aux,.(Target,GeneID), function(x) {
        if(!x$GeneID %in% tss$ensembl_gene_id) next # If gene symbol not in the TSS list
        x.tss <- tss[tss$ensembl_gene_id == x$GeneID,]
        probe <- rowRanges(getMet(data))[x$Target,]
        return(values(distanceToNearest(probe,x.tss))$distance)
      },.progress = "text",  .parallel = parallel)
      colnames(distance) <- c("Target","GeneID","distNearestTSS")
      aux <- merge(aux,distance, by = c("Target","GeneID"))
      NearGenes[[probe]] <- aux
    }
  }
  return(NearGenes)
}



