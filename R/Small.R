## Construct Multi Assay Experiment
#' @description 
#' This function will receive a gene expression and DNA methylation data objects 
#' and create a Multi Assay Experiment.
#' @param met A Summaerized Experiment, a matrix or path of rda file only containing the data.
#' @param exp A Summaerized Experiment, a matrix or path of rda file only containing the data. Rownames should be 
#' either Ensembl gene id (ensembl_gene_id) or gene symbol (external_gene_name)
#' @param genome Which is the default genome to make gene information. Options hg19 and hg38
#' @param pData A DataFrame or data.frame of the phenotype data for all participants
#' @param TCGA A logical. FALSE indicate data is not from TCGA (FALSE is default). 
#' TRUE indicates data is from TCGA and sample section will automatically filled in.
#' @return A MultiAssayExperiment object
#' @export 
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @examples
#' # NON TCGA example
#' gene.exp <- DataFrame(sample1 = c("TP53"=2.3,"PTEN"=5.4),
#'                        sample2 = c("TP53"=1.6,"PTEN"=2.3)
#' )
#' dna.met <- DataFrame(sample1 = c("cg14324200"=0.5,"cg23867494"=0.1),
#'                        sample2 =  c("cg14324200"=0.3,"cg23867494"=0.9)
#' )
#' sample.info <- DataFrame(sample.type = c("Normal", "Tumor"))
#' rownames(sample.info) <- colnames(gene.exp)
#' mae <- createMultiAssayExperiment(exp = gene.exp, met = dna.met, pData = sample.info, genome = "hg38") 
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
#'    
#'    # Consisering it is TCGA and SE
#'    mae.hg19 <- createMultiAssayExperiment(exp = exp.hg19, met =  met, TCGA = TRUE, genome = "hg19")
#'    values(getExp(mae.hg19))
#'    
#'    mae.hg38 <- createMultiAssayExperiment(exp = exp.hg38, met = met, TCGA = TRUE, genome = "hg38")
#'    values(getExp(mae.hg38))
#'    
#'    # Consisering it is TCGA and not SE
#'    mae.hg19.test <- createMultiAssayExperiment(exp = assay(exp.hg19), met =  assay(met), TCGA = TRUE, genome = "hg19")
#'    
#'    mae.hg38 <- createMultiAssayExperiment(exp = assay(exp.hg38), met = assay(met), TCGA = TRUE, genome = "hg38")
#'    values(getExp(mae.hg38))
#'    
#'    # Consisering it is not TCGA and SE
#'    # DNA methylation and gene expression Objects should have same sample names in columns
#'    not.tcga.exp <- exp.hg19 
#'    colnames(not.tcga.exp) <- substr(colnames(not.tcga.exp),1,15)
#'    not.tcga.met <- exp.hg19 
#'    colnames(not.tcga.met) <- substr(colnames(not.tcga.met),1,15)
#'    
#'    phenotype.data <- data.frame(row.names = colnames(not.tcga.exp), samples = colnames(not.tcga.exp), group = c(rep("group1",4),rep("group2",4)))
#'    mae.hg19 <- createMultiAssayExperiment(exp = not.tcga.exp, met =  not.tcga.met, TCGA = FALSE, genome = "hg19", pData = phenotype.data)
#' }
#' createMultiAssayExperiment
createMultiAssayExperiment <- function (exp, 
                                        met, 
                                        pData, 
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
  
  if(class(met) != class(as(SummarizedExperiment(),"RangedSummarizedExperiment"))){
    met <- makeSummarizedExperimentFromDNAMethylation(met, genome)
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
        message("We will consider your gene expression row names are Gene Symbols")
        extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$external_gene_name),required.cols])
        colnames(extra) <- required.cols
        values(exp) <- cbind(values(exp),extra)
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
                                pData = pData,   sampleMap = sampleMap,
                                metadata = list(TCGA= TRUE))
    
  } else {
    if(missing(pData)) stop("Please set pData argument. A data frame with samples information. All rownames should be colnames of DNA methylation and gene expression")
    # Check that we have the same number of samples
    ID <- intersect(colnames(met), colnames(exp))
    met <- met[,match(ID,colnames(met))]
    exp <- exp[,match(ID,colnames(exp))]
    pData <- pData[match(ID,rownames(pData)),,drop = FALSE]
    sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), rep("Gene expression", length(colnames(exp)))),
                           primary = c(colnames(met),colnames(exp)),
                           colname=c(colnames(met),colnames(exp)))
    if(!all(colnames(exp) == colnames(met))) stop("Please, be sure your DNA methylation matrix and gene expression matrix have the samples in the same order")
    mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                 "Gene expression" = exp),
                                pData = pData,
                                sampleMap = sampleMap,
                                metadata = list(TCGA=FALSE))
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
    exp <- makeSummarizedExperimentFromDataFrame(aux[,!grepl("external_gene_name|ensembl_gene_id",colnames(aux))],    
                                                 start.field="start_position",
                                                 end.field=c("end_position"))
    extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$ensembl_gene_id),required.cols])
    colnames(extra) <- required.cols
    values(exp) <- cbind(values(exp),extra)
  } else {
    message("We will consider your gene expression row names are Gene Symbols")
    exp$external_gene_name <- rownames(exp)
    aux <- merge(exp, gene.info, by = "external_gene_name", sort = FALSE)
    aux <- aux[!duplicated(aux$external_gene_name),]
    rownames(aux) <- aux$external_gene_name
    exp <- makeSummarizedExperimentFromDataFrame(aux[,!grepl("external_gene|ensembl_gene_id",colnames(aux))],    
                                                 start.field="start_position",
                                                 end.field=c("end_position"))
    extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$external_gene_name),required.cols])
    colnames(extra) <- required.cols
    values(exp) <- cbind(values(exp),extra)
  }
  return(exp)
}

#' @importFrom downloader download
makeSummarizedExperimentFromDNAMethylation <- function(met, genome) {
  message("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  message("Creating a SummarizedExperiment from DNA methylation input")
  if(nrow(met) > 800000) {
    plat <- "EPIC"
    annotation <- "http://zwdzwd.io/InfiniumAnnotation/current/EPIC/EPIC.manifest.rda"
    message(paste0("EPIC platform identified.\nAdding annotation from: ",annotation))
  } else {
    plat <- "450K"
    annotation <- "http://zwdzwd.io/InfiniumAnnotation/current/hm450/hm450.manifest.rda"
    message(paste0("450K platform identified.\nAdding annotation from: ",annotation))
  }
  
  if(genome == "hg38") annotation <- gsub(".rda",".hg38.rda", annotation)
  
  if(Sys.info()["sysname"] == "Windows") mode <- "wb" else  mode <- "w"
  if(!file.exists(basename(annotation))) download(annotation, basename(annotation), mode = mode)
  annotation <- get(load(basename(annotation)))
  rowRanges <- annotation[rownames(met),]
  colData <-  DataFrame(samples = colnames(met))
  assay <- data.matrix(met)
  met <- SummarizedExperiment(assays=assay,
                              rowRanges=rowRanges,
                              colData=colData)
  return(met)
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
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMultiAssayExperiment function}}.
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
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMultiAssayExperiment function}}.
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
txs <- function(genome.build = "hg19",TSS=list(upstream=NULL, downstream=NULL)){
  if(genome.build == "hg38") TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
  gene <- transcripts(Homo.sapiens, columns=c('GENEID','SYMBOL','ENSEMBL',"ENTREZID"))
  gene$GENEID <- unlist(gene$GENEID)
  gene$TXNAME <- unlist(gene$TXNAME)
  gene$SYMBOL <- unlist(gene$SYMBOL)
  gene <- gene[!is.na(gene$GENEID)]
  if(!is.null(TSS$upstream) & !is.null(TSS$downstream)) 
    gene <- promoters(gene, upstream = TSS$upstream, downstream = TSS$downstream)
  return(gene)
}

# list of TF from this paper:  http://www.sciencedirect.com/science/article/pii/S0092867410000796 
getTF <- function(genome.build = "hg19"){
  newenv <- new.env()
  data("human.TF",package = "ELMER.data",envir=newenv)
  human.TF <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  
  if(genome.build == "hg38") {
    gene <- get.GRCh("hg38", human.TF$GeneID)
  } else {
    gene <- get.GRCh("hg19",human.TF$GeneID)
  }
  gene  <- gene[!duplicated(gene),]

  return(gene)
}

#' @importFrom biomaRt getBM useMart listDatasets
get.GRCh <- function(genome="hg19", genes) {
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
  return(gene.location)
}
