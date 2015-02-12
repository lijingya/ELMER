#---------------return sample type function-------------------- -------
# Returns type: tumor, bloodNormal, tissueNormal
#according to the sampletypes
# Returns type: tumor, bloodNormal, tissueNormal
#' tcgaSampleType.
#' @param x A TCGA sample barcode.
#' @return Tissue type: Tumor, Normal, Control, TRB, cell_line, XP, XCL
#' @export 
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

#' tcgaSampleType.
#' @param tcgaId A vector list TCGA sample barcode.
#' @return standardized TCGA ID
#' @export 
standardizeTcgaId<-function(tcgaId){
  out<-substring(tcgaId,1,15)
  return(out)
}








##############-----------fetch data
## Construct MEE.data from load local data.
#' tcgaSampleType.
#' @param meth A matrix or path of rda file only containing a matrix of DNA methylation data.
#' @param exp A matrix or path of rda file only containing a matrix of expression data.
#' @param sample A data frame or path of rda file only containing sample information in data frame format.
#' @param probeInfo A GRnage object or path of rda file only containing a GRange of probe information 
#' @param geneInfo A GRnage object or path of rda file only containing a GRange of gene information (Coordinates, GENEID and SYMBOL) 
#' @param probes A vector lists probes' name. If probes are specified, the methylation and probeInfo will only contain this list of probes.
#' @param gene A vector lists genes' ID. If gene are specified, the methylation and probeInfo will only contain this list of probes.
#' @param TCGA A logical. FALSE indicate data is not from TCGA (FALSE is default). TRUE indicates data is from TCGA and sample section will automatically filled in.
#' @return MEE.data object
#' @export 
fetch.data <- function(meth,exp,sample,probeInfo,geneInfo,probes=NULL,genes=NULL,TCGA=FALSE){
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
    }
  }else if(TCGA){
    if(!is.null(meth) & is.null(exp)){
      sample <- data.frame(ID=standardizeTcgaId(colnames(meth)),meth.ID=colnames(meth),TN=sapply(colnames(meth),tcgaSampleType),stringsAsFactors = F)
    }else if(is.null(meth) & !is.null(exp)){
      sample <- data.frame(ID=standardizeTcgaId(colnames(exp)),exp.ID=colnames(exp),TN=sapply(colnames(exp),tcgaSampleType),stringsAsFactors = F)
    }else if(!is.null(meth) & !is.null(meth)){
      ID <- intersect(standardizeTcgaId(colnames(meth)),standardizeTcgaId(colnames(exp)))
      meth.ID <- colnames(meth)
      exp.ID <- colnames(exp)
      sample <- data.frame(ID=ID,meth.ID=meth.ID[match(ID,standardizeTcgaId(meth.ID))],
                           exp.ID=exp.ID[match(ID,standardizeTcgaId(exp.ID))],TN=sapply(ID,tcgaSampleType),stringsAsFactors = F)
      
    }
  }else{
    ID <- intersect(colnames(meth),colnames(exp))
    if(length(ID)==0) stop("Sample naming should consistant in methylation and expression data.")
    meth.ID <- colnames(meth)
    exp.ID <- colnames(exp)
    sample <- data.frame(ID=ID,meth.ID=meth.ID[match(ID,meth.ID)],exp.ID=exp.ID[match(ID,exp.ID)],stringsAsFactors = F)
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
    exp <- exp[sub("ID","",rownames(exp)) %in% as.character(geneInfo$GENEID),]
    geneInfo <- geneInfo[as.character(geneInfo$GENEID) %in% sub("ID","",rownames(exp))]
  } 
  if(!is.null(meth) & !is.null(exp)){
    meth <- meth[,sample$meth.ID]
    exp <- exp[,sample$exp.ID]
  }
  mee <- mee.data(meth=meth,exp=exp,sample=sample,probeInfo=probeInfo,geneInfo=geneInfo)
  return(mee)
}



