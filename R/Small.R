#---------------return sample type function-------------------- -------
# Returns type: tumor, bloodNormal, tissueNormal
#according to the sampletypes
# Returns type: tumor, bloodNormal, tissueNormal
# tcgaSampleType.
#x A TCGA sample barcode.
# Tissue type: Tumor, Normal, Control, TRB, cell_line, XP, XCL
tcgaSampleType <- function(x){
  # "Code","Definition","Short Letter Code"
  # "01","Primary solid Tumor","TP"
  # "02","Recurrent Solid Tumor","TR"
  # "03","Primary Blood Derived Cancer - Peripheral Blood","TB"
  # "04","Recurrent Blood Derived Cancer - Bone Marrow","TRBM"
  # "05","Additional - New Primary","TAP"
  # "06","Metastatic","TM"
  # "07","Additional Metastatic","TAM"
  # "08","Human Tumor Original Cells","THOC"
  # "09","Primary Blood Derived Cancer - Bone Marrow","TBM"
  # "10","Blood Derived Normal","NB"
  # "11","Solid Tissue Normal","NT"
  # "12","Buccal Cell Normal","NBC"
  # "13","EBV Immortalized Normal","NEBV"
  # "14","Bone Marrow Normal","NBM"
  # "20","Control Analyte","CELLC"
  # "40","Recurrent Blood Derived Cancer - Peripheral Blood","TRB"
  # "50","Cell Lines","CELL"
  # "60","Primary Xenograft Tissue","XP"
  # "61","Cell Line Derived Xenograft Tissue","XCL"  
  code<-substring(x,14,15)
  if (code < 10){
    out <- "Tumor"
  }else if ((code == 10) || (code == 11) || (code == 13) || (code == 14)){
    out <- "Normal"
  }else if (code==20){
    out <- "Control"
  }else if (code==40){
    out <- "TRB"
  }else if (code==50){
    out <- "cell_line"
  }else if (code==60){
    out <- "XP"
  }else if (code==61){
    out <- "XCL"
  }else{
    out <- "other"
  }
  return(out)
}

# tcgaSampleType.
# param tcgaId A vector list TCGA sample barcode.
# return Standardized TCGA IDs
standardizeTcgaId<-function(tcgaId){
  out<-substring(tcgaId,1,15)
  return(out)
}








##############-----------fetch data
## Construct MEE.data from load local data.
#' fetch.mee
#' @param meth A matrix or path of rda file only containing a matrix of 
#' DNA methylation data.
#' @param exp A matrix or path of rda file only containing a matrix of 
#' expression data.
#' @param sample A data frame or path of rda file only containing sample 
#' information in data frame format.
#' @param probeInfo A GRnage object or path of rda file only containing 
#' a GRange of probe information 
#' @param geneInfo A GRnage object or path of rda file only containing 
#' a GRange of gene information (Coordinates, GENEID and SYMBOL) 
#' @param probes A vector lists probes' name. If probes are specified, 
#' the methylation and probeInfo will only contain this list of probes.
#' @param genes A vector lists genes' ID. If gene are specified, 
#' the methylation and probeInfo will only contain this list of probes.
#' @param TCGA A logical. FALSE indicate data is not from TCGA (FALSE is default). 
#' TRUE indicates data is from TCGA and sample section will automatically filled in.
#' @return A MEE.data object
#' @export 
#' @examples
#' meth <- matrix(data=c(1:20),ncol=5,dimnames=list(paste0("probe",1:4),paste0("sample",1:5)))
#' exp <- matrix(data=c(101:110),ncol=5,dimnames=list(c("gene1","gene2"),paste0("sample",1:5)))
#' mee <- fetch.mee(meth=meth, exp=exp)
#' ## only fetch probe 1 and 3
#' mee <- fetch.mee(meth=meth, exp=exp, probes=c("probe1","probe3")) 
#' ## only fetch gene 1
#' mee <- fetch.mee(meth=meth, exp=exp, genes="gene1")
fetch.mee <- function(meth,exp,sample,probeInfo,geneInfo,probes=NULL,
                      genes=NULL,TCGA=FALSE){
  if(!missing(meth)){
    if(is.character(meth)){
      newenv <- new.env()
      load(meth, envir=newenv)
      meth <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
    }
  }else{
    meth <- NULL
  }
  
  if(!missing(exp)){
    if(is.character(exp)){
      newenv <- new.env()
      load(exp, envir=newenv)
      exp <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
    }
  }else{
    exp <- NULL
  }
  if(!missing(sample)){
    if(is.character(sample)){
      newenv <- new.env()
      load(sample, envir=newenv)
      sample <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
      sample$ID <- as.character(sample$ID)
    }
  }else if(TCGA){
    if(!is.null(meth) & is.null(exp)){
      TN <- sapply(colnames(meth),tcgaSampleType)
      TN[TN %in% "Tumor"] <- "Experiment"
      TN[TN %in% "Normal"] <- "Control"
      sample <- data.frame(ID=standardizeTcgaId(colnames(meth)),
                           meth.ID=colnames(meth),
                           TN=TN,
                           stringsAsFactors = FALSE)
    }else if(is.null(meth) & !is.null(exp)){
      TN <- sapply(colnames(exp),tcgaSampleType)
      TN[TN %in% "Tumor"] <- "Experiment"
      TN[TN %in% "Normal"] <- "Control"
      sample <- data.frame(ID=standardizeTcgaId(colnames(exp)),
                           exp.ID=colnames(exp),
                           TN=TN,
                           stringsAsFactors = FALSE)
    }else if(!is.null(meth) & !is.null(exp)){
      ID <- intersect(standardizeTcgaId(colnames(meth)),standardizeTcgaId(colnames(exp)))
      TN <- sapply(ID,tcgaSampleType)
      TN[TN %in% "Tumor"] <- "Experiment"
      TN[TN %in% "Normal"] <- "Control"
      meth.ID <- colnames(meth)
      exp.ID <- colnames(exp)
      sample <- data.frame(ID=ID,meth.ID=meth.ID[match(ID,standardizeTcgaId(meth.ID))],
                           exp.ID=exp.ID[match(ID,standardizeTcgaId(exp.ID))],
                           TN=TN,stringsAsFactors = FALSE)
      
    }
  }else if(!is.null(meth) & !is.null(exp)){
    ID <- intersect(colnames(meth),colnames(exp))
    if(length(ID)==0) 
      stop("Sample naming should consistant in methylation and expression data.")
    meth.ID <- colnames(meth)
    exp.ID <- colnames(exp)
    sample <- data.frame(ID=ID,meth.ID=meth.ID[match(ID,meth.ID)],
                         exp.ID=exp.ID[match(ID,exp.ID)],stringsAsFactors = FALSE)
  }else{
    sample <- NULL
  }
  if(!missing(probeInfo)){
    if(is.character(probeInfo)){
      newenv <- new.env()
      load(probeInfo, envir=newenv)
      probeInfo <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
    }
    if(!is.null(probes)){
      probeInfo <- probeInfo[probeInfo$name %in% probes]
    }
  }else{
    probeInfo <- NULL
  }
  
  if(!missing(geneInfo)){
    if(is.character(geneInfo)){
      newenv <- new.env()
      load(geneInfo, envir=newenv)
      geneInfo <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
    }
    if(!is.null(genes)){
      geneInfo <- geneInfo[geneInfo$GENEID %in% genes]
    }
  }else{
    geneInfo <- NULL
  }
  
  if(!is.null(probeInfo) & !is.null(meth)){
    meth <- meth[rownames(meth) %in% as.character(probeInfo$name),]
    probeInfo <- probeInfo[as.character(probeInfo$name) %in% rownames(meth)]
  } 
  if(!is.null(geneInfo) & !is.null(exp)){
    exp <- exp[rownames(exp)%in% as.character(geneInfo$GENEID),]
    geneInfo <- geneInfo[as.character(geneInfo$GENEID) %in% rownames(exp)]
  } 
  if(!is.null(meth) & !is.null(exp)){
      if(TCGA){
          meth <- meth[,sample$meth.ID]
          exp <- exp[,sample$exp.ID]
      }else{
          meth <- meth[,sample$ID]
          exp <- exp[,sample$ID]
      }
    
  }
  mee <- mee.data(meth=meth,exp=exp,sample=sample,
                  probeInfo=probeInfo,geneInfo=geneInfo)
  invisible(gc())
  return(mee)
}


#' fetch.pair
#' @param pair A data.frame or path of csv file containing pair information.
#' @param probeInfo A GRnage object or path of rda file only containing 
#' a GRange of probe information 
#' @param geneInfo A GRnage object or path of rda file only containing 
#' a GRange of gene information (Coordinates, GENEID and SYMBOL) 
#' @return A pair.data object
#' @export 
#' @examples
#' df <- data.frame(Probe=c("cg19403323","cg12213388","cg26607897"),
#' GeneID =c("ID255928","ID84451","ID55811"),
#' Symbol =c("SYT14","KIAA1804","ADCY10"),
#' Pe=c(0.003322259,0.003322259,0.003322259))
#' geneInfo <- txs()
#' ## input can be a path
#' pair <- fetch.pair(pair = df, geneInfo=geneInfo)

fetch.pair <- function(pair,probeInfo,geneInfo){
  if(!missing(pair)){
    if(is.character(pair)){
      pair <- read.csv(pair, stringsAsFactors=FALSE)
    }
  }else{
    pair <- NULL
  }
  if(!missing(probeInfo)){
    if(is.character(probeInfo)){
      newenv <- new.env()
      load(probeInfo, envir=newenv)
      probeInfo <- get(ls(newenv)[1],envir=newenv)
      # The data is in the one and only variable
    }
  }else{
    probeInfo <- NULL
  }
  
  if(!missing(geneInfo)){
    if(is.character(geneInfo)){
      newenv <- new.env()
      load(geneInfo, envir=newenv)
      geneInfo <- get(ls(newenv)[1],envir=newenv) 
      # The data is in the one and only variable
    }
  }else{
    geneInfo <- NULL
  }
  
  pair <- pair.data(pairInfo=pair,probeInfo=probeInfo,geneInfo=geneInfo)
  return(pair)
}

# splitmatix 
# @param x A matrix 
# @param by A character specify if split the matix by row or column.
# @return A list each of which is the value of each row/column in the matrix.
splitmatrix <- function(x,by="row") {
  if(by %in% "row"){
    out <- split(x, rownames(x))
  }else if (by %in% "col"){
    out <- split(x, colnames(x))
  }
  return(out)
}

#'getSymbol
#'@param mee A MEE.data or Pair object.
#'@param geneID A character which is the geneID
#'@return The gene symbol 
#'@export
#' @examples
#' geneInfo <- txs()
#' ## input can be a path
#' pair <- fetch.pair(geneInfo=geneInfo)
#' getSymbol(pair, geneID="84451")
getSymbol <- function(mee,geneID){
  gene <- unique(values(getGeneInfo(mee,geneID=geneID))[,c("GENEID","SYMBOL")])
  gene <- gene[match(geneID,gene$GENEID),"SYMBOL"]
  return(gene)
}

#'getGeneID
#'@import S4Vectors
#'@param mee A MEE.data or Pair object.
#'@param symbol A character which is the geneID
#'@return The gene ID
#'@export
#'@examples
#' geneInfo <- txs()
#' ## input can be a path
#' pair <- fetch.pair(geneInfo=geneInfo)
#' getGeneID(pair, symbol="KIAA1804")
getGeneID <- function(mee,symbol){
  gene <- unique(values(getGeneInfo(mee,symbol=symbol))[,c("GENEID","SYMBOL")])
  gene <- gene[match(symbol,gene$SYMBOL),"GENEID"]
  return(gene)
}

# binary data
# @param x A matrix.
# @param Break A value to binarize the data.
# @param Break2 A value to cut value to 3 categories.
# @return A binarized matrix.
Binary <- function(x,Break=0.3,Break2=NULL){
  if(!is.numeric(x)) stop("x need to be numeric") 
  change <- x
  if(is.null(Break2)){
    change[x > Break] <- 1
    change[x < Break | x== Break] <- 0
  }else{
    change[x < Break | x== Break] <- 0
    change[x> Break & x < Break2] <- NA
    change[x > Break2 | x== Break2] <-1 
  }
  
  return(change)    
}



# lable linear regression formula 
# @param df A data.frame object contains two variables: dependent 
# variable (Dep) and explanation variable (Exp).
# @param Dep A character specify dependent variable. The first column 
# will be dependent variable as default.
# @param Exp A character specify explanation variable. The second column 
# will be explanation variable as default.
# @return A linear regression formula
lm_eqn = function(df,Dep,Exp){
  if(missing(Dep)) Dep <- colnames(df)[1]
  if(missing(Exp)) Exp <- colnames(df)[2]
  m = lm(df[,Dep] ~ df[,Exp]);
  eq <- substitute(italic(y) == a + (b) %.% italic(x)*"\n"~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


## txs
#' txs Fetch USCS gene annotation (transcripts level) from Bioconductor package 
#' Homo.sapians. If upstream and downstream are specified in TSS list, promoter 
#' regions of USCS gene will be generated.
#' @param TSS A list contains upstream and downstream like TSS=list(upstream, downstream).
#' When upstream and downstream is specified, coordinates of promoter regions with 
#' gene annotation will be generated. 
#' @return UCSC gene annotation if TSS is not specified. Coordinates
#' of UCSC gene promoter regions with each gene annotation if TSS
#' is specified.
#' @importFrom GenomicFeatures transcripts
#' @examples
#' # get UCSC gene annotation (transcripts level)
#' \dontrun{
#' txs <- txs()
#' }
#' # get coordinate of all UCSC promoter regions
#' \dontrun{
#' txs <- txs(TSS=list(upstream=1000, downstream=1000))
#' }
#' @export
#' @import BiocGenerics GenomeInfoDb GenomicRanges
txs <- function(TSS=list(upstream=NULL, downstream=NULL)){
  gene <- transcripts(Homo.sapiens, columns=c('TXNAME','GENEID','SYMBOL'))
  gene$GENEID <- unlist(gene$GENEID)
  gene$TXNAME <- unlist(gene$TXNAME)
  gene$SYMBOL <- unlist(gene$SYMBOL)
  gene <- gene[!is.na(gene$GENEID)]
  if(!is.null(TSS$upstream) & !is.null(TSS$downstream)) 
    gene <- promoters(gene, upstream = TSS$upstream, downstream = TSS$downstream)
  return(gene)
}
