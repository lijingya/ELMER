#' @title  Construct a Multi Assay Experiment for ELMER analysis
#' @description 
#' This function will receive a gene expression and DNA methylation data objects 
#' and create a Multi Assay Experiment.
#' @param met A Summaerized Experiment, a matrix or path of rda file only containing the data.
#' @param exp A Summaerized Experiment, a matrix or path of rda file only containing the data. Rownames should be 
#' either Ensembl gene id (ensembl_gene_id) or gene symbol (external_gene_name)
#' @param genome Which is the default genome to make gene information. Options hg19 and hg38
#' @param pData A DataFrame or data.frame of the phenotype data for all participants
#' @param sampleMap  A DataFrame or data.frame of the matching samples and colnames
#'  of the gene expression and DNA methylation matrix. This should be used if your matrix
#'  have different columns names. 
#'  This object must have columns primary (sample ID) and colname (names of the columns of the matrix).
#' @param linearize.exp Take log2(exp + 1) in order to linearize relation between methylation and expression  
#' @param TCGA A logical. FALSE indicate data is not from TCGA (FALSE is default). 
#' TRUE indicates data is from TCGA and sample section will automatically filled in.
#' @param filter.probes A GRanges object contains the coordinate of probes which locate 
#'  within promoter regions or distal feature regions such as union enhancer from REMC and FANTOM5.
#'  See \code{\link{get.feature.probe}} function.
#' @param filter.genes List of genes ensemble ids to filter from object  
#' @return A MultiAssayExperiment object
#' @export 
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment makeSummarizedExperimentFromDataFrame assay assay<-
#' @examples
#' # NON TCGA example: matrices has diffetrent column names
#' gene.exp <- DataFrame(sample1 = c("TP53"=2.3,"PTEN"=5.4),
#'                        sample2 = c("TP53"=1.6,"PTEN"=2.3)
#' )
#' dna.met <- DataFrame(sample1 = c("cg14324200"=0.5,"cg23867494"=0.1),
#'                        sample2 =  c("cg14324200"=0.3,"cg23867494"=0.9)
#' )
#' sample.info <- DataFrame(sample.type = c("Normal", "Tumor"))
#' rownames(sample.info) <- colnames(gene.exp)
#' mae <- createMAE(exp = gene.exp, met = dna.met, pData = sample.info, genome = "hg38") 
#' 
#' # NON TCGA example: matrices has diffetrent column names
#' gene.exp <- DataFrame(sample1.exp = c("TP53"=2.3,"PTEN"=5.4),
#'                        sample2.exp = c("TP53"=1.6,"PTEN"=2.3)
#' )
#' dna.met <- DataFrame(sample1.met = c("cg14324200"=0.5,"cg23867494"=0.1),
#'                        sample2.met =  c("cg14324200"=0.3,"cg23867494"=0.9)
#' )
#' sample.info <- DataFrame(sample.type = c("Normal", "Tumor"))
#' rownames(sample.info) <- c("sample1","sample2")
#' sampleMap <- DataFrame(primary = c("sample1","sample1","sample2","sample2"), 
#'                       colname = c("sample1.exp","sample1.met","sample2.exp","sample2.met"))
#' mae <- createMAE(exp = gene.exp, met = dna.met, sampleMap = sampleMap, pData = sample.info, genome = "hg38") 
#' \dontrun{
#'    # TCGA example using TCGAbiolinks
#'    # Testing creating MultyAssayExperiment object
#'    # Load library
#'    library(TCGAbiolinks)
#'    library(SummarizedExperiment)
#'    
#'    samples <- c("TCGA-BA-4074", "TCGA-BA-4075", "TCGA-BA-4077", "TCGA-BA-5149",
#'                 "TCGA-UF-A7JK", "TCGA-UF-A7JS", "TCGA-UF-A7JT", "TCGA-UF-A7JV")
#'    
#'    #1) Get gene expression matrix
#'    query.exp <- GDCquery(project = "TCGA-HNSC", 
#'                          data.category = "Transcriptome Profiling", 
#'                          data.type = "Gene Expression Quantification", 
#'                          workflow.type = "HTSeq - FPKM-UQ",
#'                          barcode = samples)
#'    
#'    GDCdownload(query.exp)
#'    exp.hg38 <- GDCprepare(query = query.exp)
#'    
#'    
#'    # Aligned against Hg19
#'    query.exp.hg19 <- GDCquery(project = "TCGA-HNSC", 
#'                               data.category = "Gene expression",
#'                               data.type = "Gene expression quantification",
#'                               platform = "Illumina HiSeq", 
#'                               file.type  = "normalized_results",
#'                               experimental.strategy = "RNA-Seq",
#'                               barcode = samples,
#'                               legacy = TRUE)
#'    GDCdownload(query.exp.hg19)
#'    exp.hg19 <- GDCprepare(query.exp.hg19)
#'    
#'    # Our object needs to have emsembl gene id as rownames
#'    rownames(exp.hg19) <- values(exp.hg19)$ensembl_gene_id
#'    
#'    # DNA Methylation
#'    query.met <- GDCquery(project = "TCGA-HNSC",
#'                          legacy = TRUE,
#'                          data.category = "DNA methylation",
#'                          barcode = samples,
#'                          platform = "Illumina Human Methylation 450")
#'    
#'    GDCdownload(query.met)
#'    met <- GDCprepare(query = query.met)
#'    
#'    distal.enhancer <- get.feature.probe(genome = "hg19",platform = "450k")                             
#'    
#'    # Consisering it is TCGA and SE
#'    mae.hg19 <- createMAE(exp = exp.hg19, met =  met, TCGA = TRUE, genome = "hg19",  filter.probes = distal.enhancer)
#'    values(getExp(mae.hg19))
#'    
#'    mae.hg38 <- createMAE(exp = exp.hg38, met = met, TCGA = TRUE, genome = "hg38",  filter.probes = distal.enhancer)
#'    values(getExp(mae.hg38))
#'    
#'    # Consisering it is TCGA and not SE
#'    mae.hg19.test <- createMAE(exp = assay(exp.hg19), met =  assay(met), 
#'                               TCGA = TRUE, genome = "hg19",  
#'                               filter.probes = distal.enhancer)
#'    
#'    mae.hg38 <- createMAE(exp = assay(exp.hg38), met = assay(met), 
#'                          TCGA = TRUE, genome = "hg38",  
#'                          filter.probes = distal.enhancer)
#'    values(getExp(mae.hg38))
#'    
#'    # Consisering it is not TCGA and SE
#'    # DNA methylation and gene expression Objects should have same sample names in columns
#'    not.tcga.exp <- exp.hg19 
#'    colnames(not.tcga.exp) <- substr(colnames(not.tcga.exp),1,15)
#'    not.tcga.met <- met 
#'    colnames(not.tcga.met) <- substr(colnames(not.tcga.met),1,15)
#'    
#'    phenotype.data <- data.frame(row.names = colnames(not.tcga.exp), 
#'                                 samples = colnames(not.tcga.exp), 
#'                                 group = c(rep("group1",4),rep("group2",4)))
#'    distal.enhancer <- get.feature.probe(genome = "hg19",platform = "450k")                             
#'    mae.hg19 <- createMAE(exp = not.tcga.exp, 
#'                          met =  not.tcga.met, 
#'                          TCGA = FALSE, 
#'                          filter.probes = distal.enhancer,
#'                          genome = "hg19", 
#'                          pData = phenotype.data)
#' }
#' createMAE
createMAE <- function (exp, 
                       met, 
                       pData, 
                       sampleMap,
                       linearize.exp = FALSE,
                       filter.probes = NULL,
                       filter.genes = NULL,
                       met.platform = "450k",
                       genome = NULL,
                       TCGA = FALSE) {
  
  if(missing(genome)) stop("Please specify the genome (hg38, hg19)")
  
  # Check if input are path to rda files
  if(is.character(exp)) exp <- get(load(exp))
  if(is.character(met)) met <- get(load(met))
  
  
  # Expression data must have the ensembl_gene_id (Ensemble ID) and external_gene_name (Gene Symbol)
  required.cols <- c("external_gene_name", "ensembl_gene_id")
  # If my input is a data frame we will need to add metadata information for the ELMER analysis steps
  if(class(exp) != class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
    exp <- makeSummarizedExperimentFromGeneMatrix(exp, genome)
  }
  # Add this here ?
  if(linearize.exp) assay(exp) <- log2(assay(exp) + 1)
  
  if(class(met) != class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
    met <- makeSummarizedExperimentFromDNAMethylation(met, genome, met.platform)
  }
  
  # Select the regions from DNA methylation that overlaps enhancer.
  if(!is.null(filter.probes)){
    if(is.character(filter.probes)){
      filter.probes <- get(load(filter.probes))
    }
  } 
  if(!is.null(filter.probes) & !is.null(met)){
    met <- met[rownames(met) %in% names(filter.probes),]
  }
  if(!is.null(filter.genes) & !is.null(exp)){
    exp <- exp[rownames(exp) %in% names(filter.genes),]
  } 
  
  # We will need to check if the fields that we need exists.
  # Otherwise we will need to create them
  if(class(exp) == class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
    required.cols <- required.cols[!required.cols %in% colnames(values(exp))]
    if(length(required.cols) > 0) {
      gene.info <- TCGAbiolinks:::get.GRCh.bioMart(genome)
      colnames(gene.info)[grep("external_gene", colnames(gene.info))] <- "external_gene_name"
      if(all(grepl("ENSG",rownames(exp)))) {
        extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$ensembl_gene_id),required.cols])
        colnames(extra) <- required.cols
        values(exp) <- cbind(values(exp),extra)
      } else {
        stop("Please the gene expression matrix should receive ENSEMBLE IDs")
        #extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$external_gene_name),required.cols])
        #colnames(extra) <- required.cols
        #values(exp) <- cbind(values(exp),extra)
      }
    }
  } 
  if(TCGA){
    message("Checking samples have both DNA methylation and Gene expression and they are in the same order...")
    # If it is not TCGA we will assure the sample has both DNA methylation and gene expression
    ID <- intersect(substr(colnames(met),1,15), substr(colnames(exp),1,15))
    
    # Get only samples with both DNA methylation and Gene expression
    met <- met[,match(ID,substr(colnames(met),1,15))]
    exp <- exp[,match(ID,substr(colnames(exp),1,15))]
    stopifnot(all(substr(colnames(exp),1,15) == substr(colnames(met),1,15)))
    
    # Get clinical information
    if(missing(pData)) {
      pData <- TCGAbiolinks:::colDataPrepare(c(colnames(met), colnames(exp)))
      sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), rep("Gene expression", length(colnames(exp)))),
                             primary = substr(c(colnames(met),colnames(exp)),1,16),
                             colname=c(colnames(met),colnames(exp)))
      pData$barcode <- NULL
      pData <- pData[!duplicated(pData),]      
      rownames(pData) <- pData$sample
    }
    message("Creating MultiAssayExperiment")
    mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                 "Gene expression" = exp),
                                pData = pData,   
                                sampleMap = sampleMap,
                                metadata = list(TCGA= TRUE, genome = genome))
  } else {
    
    if(missing(pData)){
      message <- paste("Please set pData argument. A data frame with samples", 
                       "information. All rownames should be colnames of DNA",
                       "methylation and gene expression. An example is showed",
                       "in MultiAssayExperiment documentation",
                       "(access it with ?MultiAssayExperiment)")
      stop(message)
    }
    
    if(missing(sampleMap)){
      # Check that we have the same number of samples
      message("Removing samples not found in both DNA methylation and gene expression (we are considering the names of the gene expression and DNA methylation columns to be the same) ")
      ID <- intersect(colnames(met), colnames(exp))
      met <- met[,match(ID,colnames(met))]
      exp <- exp[,match(ID,colnames(exp))]
      
      if(!all(colnames(exp) == colnames(met))) 
        stop("Error DNA methylation matrix and gene expression matrix are not in the same order")
      
      pData <- pData[match(ID,rownames(pData)),,drop = FALSE]
      sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), 
                                      rep("Gene expression", length(colnames(exp)))),
                             primary = c(colnames(met),colnames(exp)),
                             colname=c(colnames(met),colnames(exp)))
      mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                   "Gene expression" = exp),
                                  pData = pData,
                                  sampleMap = sampleMap,
                                  metadata = list(TCGA=FALSE, genome = genome))
    } else {
      # Check that we have the same number of samples
      if(!all(c("primary","colname") %in% colnames(sampleMap))) 
        stop("sampleMap should have the following columns: primary (sample ID) and colname(DNA methylation and gene expression sample [same as the colnames of the matrix])")
      if(!any(rownames(pData) %in% sampleMap$primary))
        stop("pData row names should be mapped to sampleMap primary column ")
      # Find which samples are DNA methylation and gene expression
      sampleMap.met <- sampleMap[sampleMap$colname %in% colnames(met),,drop = FALSE]
      sampleMap.exp <- sampleMap[sampleMap$colname %in% colnames(exp),,drop = FALSE]
      
      # Which ones have both DNA methylation and gene expresion?
      commun.samples <- intersect(sampleMap.met$primary,sampleMap.exp$primary)
      
      # Remove the one that does not have both data
      sampleMap.met <- sampleMap.met[match(sampleMap.met$primary,commun.samples),,drop = FALSE]
      sampleMap.exp <- sampleMap.exp[match(sampleMap.exp$primary,commun.samples),,drop = FALSE]
      
      # Ordering samples to be matched
      met <- met[,sampleMap.met$colname,drop = FALSE]
      exp <- exp[,sampleMap.exp$colname,drop = FALSE]
      
      if(!all(sampleMap.met$primary == sampleMap.exp$primary)) 
        stop("Error DNA methylation matrix and gene expression matrix are not in the same order")
      
      pData <- pData[match(commun.samples,rownames(pData)),,drop = FALSE]
      sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), 
                                      rep("Gene expression", length(colnames(exp)))),
                             primary = commun.samples,
                             colname=c(colnames(met),colnames(exp)))
      mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                   "Gene expression" = exp),
                                  pData = pData,
                                  sampleMap = sampleMap,
                                  metadata = list(TCGA=FALSE, genome = genome))
    }
  }
  return(mae)
}

makeSummarizedExperimentFromGeneMatrix <- function(exp, genome = genome){
  message("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  message("Creating a SummarizedExperiment from gene expression input")
  gene.info <- TCGAbiolinks:::get.GRCh.bioMart(genome)
  gene.info$chromosome_name <- paste0("chr",gene.info$chromosome_name)
  colnames(gene.info)[grep("external_gene", colnames(gene.info))] <- "external_gene_name"
  gene.info$strand[gene.info$strand == 1] <- "+"
  gene.info$strand[gene.info$strand == -1] <- "-"
  exp <- as.data.frame(exp)
  required.cols <- c("external_gene_name", "ensembl_gene_id")
  
  if(all(grepl("ENSG",rownames(exp)))) {
    exp$ensembl_gene_id <- rownames(exp)
    aux <- merge(exp, gene.info, by = "ensembl_gene_id", sort = FALSE)
    aux <- aux[!duplicated(aux$ensembl_gene_id),]
    rownames(aux) <- aux$ensembl_gene_id
    aux$entrezgene <- NULL
    exp <- makeSummarizedExperimentFromDataFrame(aux[,!grepl("external_gene_name|ensembl_gene_id",colnames(aux))],    
                                                 start.field="start_position",
                                                 end.field=c("end_position"))
    extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$ensembl_gene_id),required.cols])
    colnames(extra) <- required.cols
    values(exp) <- cbind(values(exp),extra)
  } else {
    stop("Please the gene expression matrix should receive ENSEMBLE IDs (ENSG)")
  }
  return(exp)
}

#' @importFrom downloader download
makeSummarizedExperimentFromDNAMethylation <- function(met, genome, met.platform) {
  message("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  message("Creating a SummarizedExperiment from DNA methylation input")
  
  # Instead of looking on the size, it is better to set it as a argument as the annotation is different
  annotation <-   getInfiniumAnnotation(met.platform, genome)
  rowRanges <- annotation[rownames(met),,drop=FALSE]
  
  # Remove masked probes, besed on the annotation
  rowRanges <- rowRanges[!rowRanges$MASK.mapping]
  
  colData <-  DataFrame(samples = colnames(met))
  met <- met[rownames(met) %in% names(rowRanges),,drop = FALSE]
  assay <- data.matrix(met)
  met <- SummarizedExperiment(assays=assay,
                              rowRanges=rowRanges,
                              colData=colData)
  return(met)
}

getInfiniumAnnotation <- function(plat = "450K", genome = "hg38"){
  if(plat == "EPIC") {
    annotation <- "http://zwdzwd.io/InfiniumAnnotation/current/EPIC/EPIC.manifest.rda"
  } else {
    annotation <- "http://zwdzwd.io/InfiniumAnnotation/current/hm450/hm450.manifest.rda"
  }
  if(genome == "hg38") annotation <- gsub(".rda",".hg38.rda", annotation)
  message(paste0("Adding annotation for DNA methylation from: ",annotation))
  
  if(!file.exists(basename(annotation))) {
    if(Sys.info()["sysname"] == "Windows") mode <- "wb" else  mode <- "w"
    download(annotation, basename(annotation), mode = mode)
  }
  annotation <- get(load(basename(annotation)))
  return(annotation)  
}

#' fetch.pair to generate Pair class object.
#' @description 
#' fetch.pair is a funtion to take in enhancer-gene linkage prediction information,
#' probe information and gene annotation generating a Pair class object, which is the 
#' input for plotting functions. Options (pair, probeInfo, geneInfo) can
#' take in R object or read files by specifying file paths. 
#' @param pair A data.frame (R object) or a path of XX.csv file containing pair information such as
#' output of function \code{\link{get.pair}}.
#' @param probeInfoA GRnage object or a path of XX.rda file which only contains a GRange of probe information.
#' @param geneInfo A GRnage object or path of XX.rda file which only contains a GRange of gene 
#' information such as Coordinates, GENEID and SYMBOL. 
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

#'getSymbol to report gene symbol from id
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE function}}.
#' @param geneID A character which is the ensembl_gene_id
#' @return The gene symbol for input genes.
#' @export
#' @examples
#' data(elmer.data.example)
#' getSymbol(data, geneID="ENSG00000143067")
getSymbol <- function(data,geneID){
  gene <- unique(values(getExp(data))[,c("ensembl_gene_id","external_gene_name")])
  gene <- gene[match(geneID,gene$ensembl_gene_id),"external_gene_name"]
  return(gene)
}

#'getGeneID to report gene id from symbol
#'@importFrom S4Vectors values
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE function}}.
#'@param symbol A vector of characters which are gene symbols 
#'@return The gene ID for these gene symbols
#'@export
#'@examples
#' data(elmer.data.example)
#' getGeneID(data, symbol="ZNF697")
getGeneID <- function(data,symbol){
  gene <- unique(values(getExp(data))[,c("ensembl_gene_id","external_gene_name")])
  gene <- gene[match(symbol,gene$external_gene_name),"ensembl_gene_id"]
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


#' txs to fetch USCS gene annotation (transcripts level) from Bioconductor package Homo.sapians. 
#' If upstream and downstream are specified in TSS list, promoter regions of USCS gene will be generated.
#' @description 
#' txs is a function to fetch USCS gene annotation (transcripts level) from Bioconductor package Homo.sapians.
#' If upstream and downstream are specified in TSS list, promoter regions of USCS gene will be generated.
#' @param TSS A list. Contains upstream and downstream like TSS=list(upstream, downstream).
#'  When upstream and downstream is specified, coordinates of promoter regions with gene annotation will be generated.
#' @param genome.build Use TxDb.Hsapiens.UCSC.hg38.knownGene instead of TxDb.Hsapiens.UCSC.hg19.knownGene.
#' Options: hg19 (default) and hg38.
#' @return UCSC gene annotation if TSS is not specified. Coordinates of UCSC gene promoter regions if TSS is specified.
#' @examples
#' # get UCSC gene annotation (transcripts level)
#' \dontrun{
#'     txs <- txs()
#'     txs <- txs(genome.build = "hg38")
#' }
#' # get coordinate of all UCSC promoter regions +/-1000bp of TSSs
#' \dontrun{
#' txs <- txs(TSS=list(upstream=1000, downstream=1000))
#' }
#' @export
#' @author Lijing Yao (maintainer: lijingya@usc.edu)
#' @import GenomeInfoDb
#' @importFrom GenomicFeatures transcripts
#' @importFrom rvest %>%
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene Homo.sapiens
txs <- function(genome = "hg38",TSS=list(upstream=NULL, downstream=NULL)){
  if(genome == "hg38") TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
  gene <- suppressMessages(transcripts(Homo.sapiens, columns=c('GENEID','SYMBOL','ENSEMBL',"ENTREZID")))
  gene$GENEID <- unlist(gene$GENEID)
  gene$TXNAME <- unlist(gene$TXNAME)
  gene$SYMBOL <- unlist(gene$SYMBOL)
  gene <- gene[!is.na(gene$GENEID)]
  if(!is.null(TSS$upstream) & !is.null(TSS$downstream)) 
    gene <- promoters(gene, upstream = TSS$upstream, downstream = TSS$downstream)
  return(gene)
}

#' getTSS to fetch GENCODE gene annotation (transcripts level) from Bioconductor package biomaRt
#' If upstream and downstream are specified in TSS list, promoter regions of GENCODE gene will be generated.
#' @description 
#' getTSS to fetch GENCODE gene annotation (transcripts level) from Bioconductor package biomaRt
#' If upstream and downstream are specified in TSS list, promoter regions of GENCODE gene will be generated.
#' @param TSS A list. Contains upstream and downstream like TSS=list(upstream, downstream).
#'  When upstream and downstream is specified, coordinates of promoter regions with gene annotation will be generated.
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @return GENCODE gene annotation if TSS is not specified. Coordinates of GENCODE gene promoter regions if TSS is specified.
#' @examples
#' # get UCSC gene annotation (transcripts level)
#' \dontrun{
#'     txs <- txs()
#'     txs <- txs(genome.build = "hg38")
#' }
#' # get coordinate of all UCSC promoter regions +/-1000bp of TSSs
#' \dontrun{
#' txs <- txs(TSS=list(upstream=1000, downstream=1000))
#' }
#' @export
#' @author Lijing Yao (maintainer: lijingya@usc.edu)
#' @import GenomeInfoDb
#' @importFrom GenomicFeatures transcripts
#' @importFrom rvest %>%
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene Homo.sapiens
getTSS <- function(genome="hg38",TSS=list(upstream=NULL, downstream=NULL)){
  if (genome == "hg19"){
    # for hg19
    ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                       host = "feb2014.archive.ensembl.org",
                       path = "/biomart/martservice" ,
                       dataset = "hsapiens_gene_ensembl")
    attributes <- c("chromosome_name",
                    "start_position",
                    "end_position", "strand",
                    "ensembl_transcript_id",
                    "ensembl_gene_id", "entrezgene",
                    "external_gene_id")
  } else {
    # for hg38
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    attributes <- c("chromosome_name",
                    "start_position",
                    "end_position", "strand",
                    "ensembl_gene_id", 
                    "ensembl_transcript_id",
                    "external_gene_name")
  }
  chrom <- c(1:22, "X", "Y","M","*")
  description <- listDatasets(ensembl)[listDatasets(ensembl)$dataset=="hsapiens_gene_ensembl",]$description
  message(paste0("Downloading genome information. Using: ", description))
  
  filename <-  paste0(gsub("[[:punct:]]| ", "_",description),"_tss.rda")
  if(!file.exists(filename)) {
    tss <- getBM(attributes = attributes, filters = c("chromosome_name"), values = list(chrom), mart = ensembl)
    tss <- tss[!duplicated(tss$ensembl_transcript_id),]
    save(tss, file = filename)
  } else {
    tss <- get(load(filename))  
  } 
  tss$chromosome_name <-  paste0("chr", tss$chromosome_name)
  tss$strand[tss$strand == 1] <- "+" 
  tss$strand[tss$strand == -1] <- "-" 
  tss <- makeGRangesFromDataFrame(tss,start.field = "start_position", end.field = "end_position", keep.extra.columns = TRUE)
  if(genome == "hg19") tss$external_gene_name <- tss$external_gene_id
  if(!is.null(TSS$upstream) & !is.null(TSS$downstream)) 
    tss <- promoters(tss, upstream = TSS$upstream, downstream = TSS$downstream)
  
  return(tss)
}

#' @title  Get human TF list from the UNiprot database
#' @description This function gets the last version of human TF list from the UNiprot database
#' @importFrom readr read_tsv
#' @return A data frame with the ensemble gene id and entrezgene and gene symbol.
getTF <- function(genome.build = "hg38"){
  uniprotURL <- "http://www.uniprot.org/uniprot/?"
  query <- "query=reviewed:yes+AND+organism:9606+AND+%22transcription+factor%22&sort=score"
  fields <- "columns=id,entry%20name,protein%20names,genes,database(GeneWiki),database(Ensembl),database(GeneID)"
  format <- "format=tab"
  human.TF <- readr::read_tsv(paste0(uniprotURL,
                                     paste(query, fields,format, sep = "&")),
                              col_types = "ccccccc") 
  gene <- get.GRCh(genome.build, gsub(";","",human.TF$`Cross-reference (GeneID)`))
  gene  <- gene[!duplicated(gene),]
  return(gene)
}

#' @importFrom biomaRt getBM useMart listDatasets
get.GRCh <- function(genome="hg38", genes) {
  if (genome == "hg19"){
    # for hg19
    ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                       host = "feb2014.archive.ensembl.org",
                       path = "/biomart/martservice" ,
                       dataset = "hsapiens_gene_ensembl")
    attributes <- c("ensembl_gene_id", "entrezgene","external_gene_id")
  } else {
    # for hg38
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    attributes <- c("ensembl_gene_id", "entrezgene","external_gene_name")
  }
  gene.location <- getBM(attributes = attributes,
                         filters = c("entrezgene"),
                         values = list(genes), mart = ensembl)
  colnames(gene.location) <-  c("ensembl_gene_id", "entrezgene","external_gene_name")
  gene.location <- gene.location[match(genes,gene.location$entrezgene),]
  return(gene.location)
}


# To find for each probe the know motif we will use HOMER software (http://homer.salk.edu/homer/)
# Step:
# 1 - get DNA methylation probes annotation with the regions
# 2 - Make a bed file from it
# 3 - Execute section: Finding Instance of Specific Motifs from http://homer.salk.edu/homer/ngs/peakMotifs.html to the HOCOMOCO TF motifs
# Also, As HOMER is using more RAM than the available we will split the files in to 100k probes.
# Obs: for each probe we create a winddow of 500 bp (-size 500) around it. This might lead to false positives, but will not have false negatives.
# The false posives will be removed latter with some statistical tests.
getBedForDNAmethylation <- function(){
  TFBS.motif <- "http://hocomoco.autosome.ru/final_bundle/HUMAN/mono/HOCOMOCOv10_HUMAN_mono_homer_format_0.0001.motif"
  if(!file.exists(basename(TFBS.motif))) downloader::download(TFBS.motif,basename(TFBS.motif))
  for(plat in c("450K","EPIC")){
    for(gen in c("hg19","hg38")){
      
      file <- paste0(plat,gen,".txt")
      print(file)
      if(!file.exists(file)){
        # STEP 1
        gr <- getInfiniumAnnotation(plat = plat,genome =  gen)
        
        # This will remove masked probes. They have poor quality and might be arbitrarily positioned (Wanding Zhou)
        print(table(gr$MASK.mapping))
        gr <- gr[!gr$MASK.mapping]
        print(table(gr$MASK.mapping))
        
        df <- data.frame(seqnames=seqnames(gr),
                         starts=as.integer(start(gr)),
                         ends=end(gr),
                         names=names(gr),
                         scores=c(rep(".", length(gr))),
                         strands=strand(gr))
        step <- 100000 # nb of lines in each file 100K was selected to not explode RAM
        n <- nrow(df)
        for(j in 0:floor(n/step)){
          # STEP 2
          file.aux <- paste0(plat,gen,"_",j,".bed")
          if(!file.exists(gsub(".bed",".txt",file.aux))){
            end <- ifelse(((j + 1) * step) > n, n,((j + 1) * step))
            write.table(df[((j * step) + 1):end,], file = file.aux, col.names = F, quote = F,row.names = F,sep = "\t")
            
            # STEP 3
            cmd <- paste0("source ~/.bash_rc; annotatePeaks.pl " ,file.aux, " ", gen, " -m ", basename(TFBS.motif), " -size 500 -cpu 12 > ", gsub(".bed",".txt",file.aux))
            system(cmd)
          }
        }
        # We will merge the results from each file into one
        peaks <- NULL
        for(j in 0:floor(n/step)){
          aux <-  read_tsv(paste0(plat,gen,"_",j,".txt"))
          colnames(aux)[1] <- "PeakID"
          if(is.null(peaks)) {
            peaks <- aux
          } else {
            peaks <- rbind(peaks, aux)
          }
        }
        write_tsv(peaks,path=file,col_names = TRUE)
        print("DONE!")
        gc()
      }
    }
  }
}

# This code will read the table with the motifs, save it as a sparce matrix
# and save all as a .rda that will be placed in ELMER.data
# Should this code be moved to ELMER.data?
prepare_object <- function(){
  # command 
  # annotatePeaks.pl /Users/chedraouisil/ELMER/450Khg19.bed hg19 -m HOCOMOCOv10_HUMAN_mono_homer_format_0.0001.motif > 450hg19.txt
  for(plat in c("450K","EPIC")){
    for(gen in c("hg19","hg38")){
      file <- paste0(plat,gen,".txt")
      motifs <- readr::read_tsv(file)
      # From 1 to 21 we have annotations then we have 640 motifs
      matrix <- Matrix::Matrix(0, nrow = nrow(motifs), ncol = 640,sparse = TRUE)
      print(object.size(matrix))
      colnames(matrix) <- gsub(" Distance From Peak\\(sequence,strand,conservation\\)","",colnames(motifs)[-c(1:21)])
      rownames(matrix) <- motifs$PeakID
      print(object.size(matrix))
      matrix[!is.na(motifs[,-c(1:21)])] <- 1
      matrix <- as(matrix, "nsparseMatrix")
      assign(paste0("Probes.motif.",gen,".",plat),matrix) # For each probe that there is a bind to the motif, we will add as 1 in the matrix. (Are hg38, hg19,450k,epic matrices equal?)
      rm(matrix)
      rm(motifs)
      gc()
    }
  }
  save(Probes.motif.hg19.450K, file = "Probes.motif.hg19.450K.rda", compress = "xz")
  save(Probes.motif.hg38.450K, file = "Probes.motif.hg38.450K.rda", compress = "xz")
  save(Probes.motif.hg19.EPIC, file = "Probes.motif.hg19.EPIC.rda", compress = "xz")
  save(Probes.motif.hg38.EPIC, file = "Probes.motif.hg38.EPIC.rda", compress = "xz")
  
}

#' @title Get family of transcription factors
#' @description This function will use TF Class database to create the object
#' that maps for each TF the members of its family. TF in the same family have 
#' high correlared PWM.
#' @importFrom rvest html_table
#' @importFrom xml2 read_html 
#' @return A list of TFs and its family members
createMotifRelevantTfs <- function(){
  if(!file.exists("motif.relavent.TFs.rda")){
    # Download from http://hocomoco.autosome.ru/human/mono
    tf.family <- "http://hocomoco.autosome.ru/human/mono" %>% read_html()  %>%  html_table()
    tf.family <- tf.family[[1]]
    # Split TF for each family, this will help us map for each motif which are the some ones in the family
    # basicaly: for a TF get its family then get all TF in that family
    family <- split(tf.family,f = tf.family$`TF family`)
    motif.relavent.TFs <- plyr::alply(tf.family,1, function(x){  
      f <- x$`TF family`
      if(f == "") return(x$`Transcription factor`) # Casse without family, we will get only the same object
      return(unique(family[as.character(f)][[1]]$`Transcription factor`))
    },.progress = "text")
    #names(motif.relavent.TFs) <- tf.family$`Transcription factor`
    names(motif.relavent.TFs) <- tf.family$Model
    # Cleaning object
    attr(motif.relavent.TFs,which="split_type") <- NULL
    attr(motif.relavent.TFs,which="split_labels") <- NULL
    save(motif.relavent.TFs, file = "motif.relavent.TFs.rda")
  } else {
    motif.relavent.TFs <- get(load("motif.relavent.TFs.rda"))
  }
  return(motif.relavent.TFs)
}


