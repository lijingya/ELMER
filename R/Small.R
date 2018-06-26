#' @title  Construct a Multi Assay Experiment for ELMER analysis
#' @description 
#' This function will receive a gene expression and DNA methylation data objects 
#' and create a Multi Assay Experiment.
#' @param met A Summaerized Experiment, a matrix or path of rda file only containing the data.
#' @param exp A Summaerized Experiment, a matrix or path of rda file only containing the data. Rownames should be 
#' either Ensembl gene id (ensembl_gene_id) or gene symbol (external_gene_name)
#' @param genome Which is the default genome to make gene information. Options hg19 and hg38
#' @param colData A DataFrame or data.frame of the phenotype data for all participants
#' @param sampleMap  A DataFrame or data.frame of the matching samples and colnames
#'  of the gene expression and DNA methylation matrix. This should be used if your matrix
#'  have different columns names. 
#'  This object must have columns primary (sample ID) and colname (names of the columns of the matrix).
#' @param linearize.exp Take log2(exp + 1) in order to linearize relation between methylation and expression  
#' @param met.platform DNA methylation platform "450K" or "EPIC"
#' @param TCGA A logical. FALSE indicate data is not from TCGA (FALSE is default). 
#' TRUE indicates data is from TCGA and sample section will automatically filled in.
#' @param filter.probes A GRanges object contains the coordinate of probes which locate 
#'  within promoter regions or distal feature regions such as union enhancer from REMC and FANTOM5.
#'  See \code{\link{get.feature.probe}} function.
#' @param filter.genes List of genes ensemble ids to filter from object  
#' @param save If TRUE, MAE object will be saved into a file named as the argument save.file if this was set, otherwise as mae_genome_met.platform.rda.
#' @param save.filename Name of the rda file to save the object (must end in .rda)
#' @param met.na.cut Define the percentage of NA that the line should have to remove the probes for humanmethylation platforms.
#' @return A MultiAssayExperiment object
#' @export 
#' @importFrom MultiAssayExperiment MultiAssayExperiment 
#' @importFrom SummarizedExperiment SummarizedExperiment makeSummarizedExperimentFromDataFrame assay assay<-
#' @examples
#' # NON TCGA example: matrices has diffetrent column names
#' gene.exp <- S4Vectors::DataFrame(sample1.exp = c("ENSG00000141510"=2.3,"ENSG00000171862"=5.4),
#'                   sample2.exp = c("ENSG00000141510"=1.6,"ENSG00000171862"=2.3))
#' dna.met <- S4Vectors::DataFrame(sample1.met = c("cg14324200"=0.5,"cg23867494"=0.1),
#'                        sample2.met =  c("cg14324200"=0.3,"cg23867494"=0.9))
#' sample.info <- S4Vectors::DataFrame(primary =  c("sample1","sample2"), 
#'                                     sample.type = c("Normal", "Tumor"))
#' sampleMap <- S4Vectors::DataFrame(primary = c("sample1","sample1","sample2","sample2"), 
#'                                   colname = c("sample1.exp","sample1.met",
#'                                               "sample2.exp","sample2.met"))
#' mae <- createMAE(exp = gene.exp, 
#'                  met = dna.met, 
#'                  sampleMap = sampleMap, 
#'                  met.platform ="450K",
#'                  colData = sample.info, 
#'                  genome = "hg38") 
#' # You can also use sample Mapping and Sample information tables from a tsv file
#' # You can use the createTSVTemplates function to create the tsv files
#' readr::write_tsv(as.data.frame(sampleMap), path = "sampleMap.tsv")
#' readr::write_tsv(as.data.frame(sample.info), path = "sample.info.tsv")
#' mae <- createMAE(exp = gene.exp, 
#'                  met = dna.met, 
#'                  sampleMap = "sampleMap.tsv", 
#'                  met.platform ="450K",
#'                  colData = "sample.info.tsv", 
#'                  genome = "hg38") 
#' 
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
#'    mae.hg19 <- createMAE(exp = exp.hg19, 
#'                          met =  met, 
#'                          TCGA = TRUE, 
#'                          genome = "hg19",  
#'                          filter.probes = distal.enhancer)
#'    values(getExp(mae.hg19))
#'    
#'    mae.hg38 <- createMAE(exp = exp.hg38, met = met, 
#'                         TCGA = TRUE, genome = "hg38",  
#'                         filter.probes = distal.enhancer)
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
#'                          colData = phenotype.data)
#' }
#' createMAE
createMAE <- function (exp, 
                       met, 
                       colData, 
                       sampleMap,
                       linearize.exp = FALSE,
                       filter.probes = NULL,
                       met.na.cut = 0.2,
                       filter.genes = NULL,
                       met.platform = "450K",
                       genome = NULL,
                       save = TRUE,
                       save.filename,
                       TCGA = FALSE) {
  
  if(missing(genome)) stop("Please specify the genome (hg38, hg19)")
  
  # Check if input are path to rda files
  if(is.character(exp)) exp <- get(load(exp))
  if(is.character(met)) met <- get(load(met))
  
  suppressMessages({
    
    if(!missing(colData)) { 
      if(is.character(colData)) { 
        colData <- as.data.frame(read_tsv(colData))
        rownames(colData) <- colData$primary
      }
    }
    if(!missing(sampleMap)) { 
      if(is.character(sampleMap)) sampleMap <- read_tsv(sampleMap)
    }
  })  
  
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
  met <- met[rowMeans(is.na(assay(met))) < met.na.cut, ]
  
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
      }
    }
  } 
  if(TCGA){
    message("Checking if samples have both DNA methylation and Gene expression and if they are in the same order...")
    # If it is not TCGA we will assure the sample has both DNA methylation and gene expression
    ID <- intersect(substr(colnames(met),1,16), substr(colnames(exp),1,16))
    
    # Get only samples with both DNA methylation and Gene expression
    met <- met[,match(ID,substr(colnames(met),1,16))]
    exp <- exp[,match(ID,substr(colnames(exp),1,16))]
    stopifnot(all(substr(colnames(exp),1,16) == substr(colnames(met),1,16)))
    stopifnot(ncol(exp) == ncol(met))
    
    # Get clinical information
    if(missing(colData)) {
      colData <- TCGAbiolinks:::colDataPrepare(c(colnames(met), colnames(exp)))
      # This will keep the same strategy the old ELMER version used:
      # Every type of tumor samples (starts with T) will be set to tumor and
      # every type of normal samples   (starts with N) will be set to normal 
      # See : https://github.com/lijingya/ELMER/blob/3e050462aa41c8f542530ccddc8fa607207faf88/R/Small.R#L8-L48
      colData$TN <- NA
      colData[grep("^N",colData$shortLetterCode),"TN"] <- "Normal" 
      colData[grep("^T",colData$shortLetterCode),"TN"] <- "Tumor" 
      
      colData$barcode <- NULL
      colData <- colData[!duplicated(colData),]      
      rownames(colData) <- colData$sample
    } 
    if(missing(sampleMap)) {
      sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), rep("Gene expression", length(colnames(exp)))),
                             primary = substr(c(colnames(met),colnames(exp)),1,16),
                             colname=c(colnames(met),colnames(exp)))
    }
    
    message("Creating MultiAssayExperiment")
    mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                 "Gene expression" = exp),
                                colData = colData,   
                                sampleMap = sampleMap,
                                metadata = list(TCGA= TRUE, genome = genome, met.platform = met.platform ))
  } else {
    
    if(missing(colData)){
      message <- paste("Please set colData argument. A data frame with samples", 
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
      
      colData <- colData[match(ID,rownames(colData)),,drop = FALSE]
      sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), 
                                      rep("Gene expression", length(colnames(exp)))),
                             primary = c(colnames(met),colnames(exp)),
                             colname=c(colnames(met),colnames(exp)))
      mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                   "Gene expression" = exp),
                                  colData = colData,
                                  sampleMap = sampleMap,
                                  metadata = list(TCGA=FALSE, genome = genome, met.platform = met.platform ))
    } else {
      # Check that we have the same number of samples
      if(!all(c("primary","colname") %in% colnames(sampleMap))) 
        stop("sampleMap should have the following columns: primary (sample ID) and colname(DNA methylation and gene expression sample [same as the colnames of the matrix])")
      #if(!any(rownames(colData) %in% sampleMap$primary))
      #  stop("colData row names should be mapped to sampleMap primary column ")
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
      
      colData <- colData[match(commun.samples,rownames(colData)),,drop = FALSE]
      sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), 
                                      rep("Gene expression", length(colnames(exp)))),
                             primary = commun.samples,
                             colname=c(colnames(met),colnames(exp)))
      mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                   "Gene expression" = exp),
                                  colData = colData,
                                  sampleMap = sampleMap,
                                  metadata = list(TCGA=FALSE, genome = genome, met.platform = met.platform ))
    }
  }
  if(save) {
    if(missing(save.filename)) save.filename <- paste0("mae_",genome,"_",met.platform,".rda")
    save(mae, file = save.filename,compress = "xz")
    message("MAE saved as: ", save.filename)
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
#' @importFrom S4Vectors DataFrame
makeSummarizedExperimentFromDNAMethylation <- function(met, genome, met.platform) {
  message("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  message("Creating a SummarizedExperiment from DNA methylation input")
  
  # Instead of looking on the size, it is better to set it as a argument as the annotation is different
  annotation <-   getInfiniumAnnotation(met.platform, genome)
  
  rowRanges <- annotation[names(annotation) %in% rownames(met),,drop=FALSE]
  
  # Remove masked probes, besed on the annotation
  rowRanges <- rowRanges[!rowRanges$MASK_general]
  
  colData <-  DataFrame(samples = colnames(met))
  met <- met[rownames(met) %in% names(rowRanges),,drop = FALSE]
  met <- met[names(rowRanges),,drop = FALSE]
  assay <- data.matrix(met)
  met <- SummarizedExperiment(assays=assay,
                              rowRanges=rowRanges,
                              colData=colData)
  return(met)
}

getInfiniumAnnotation <- function(plat = "450K", genome = "hg38"){
  if(tolower(genome) == "hg19" & toupper(plat) == "450K" ) return(getdata("hm450.hg19.manifest"))
  if(tolower(genome) == "hg19" & toupper(plat) == "EPIC" ) return(getdata("EPIC.hg19.manifest"))
  if(tolower(genome) == "hg38" & toupper(plat) == "450K" ) return(getdata("hm450.hg38.manifest"))
  if(tolower(genome) == "hg38" & toupper(plat) == "EPIC" ) return(getdata("EPIC.hg38.manifest"))
}

getdata <- function(...)
{
  e <- new.env()
  name <- data(..., package = "ELMER.data",envir = e)[1]
  e[[ls(envir = e)[1]]]
}

#' Create examples files for Sample mapping and information used in createMAE function
#' @description 
#' This function will receive the DNA methylation and gene expression matrix and will create
#' some examples of table for the argument colData and sampleMap used in ceeateMae function.
#' @param met DNA methylation matrix or Summarized Experiment
#' @param exp Gene expression matrix or Summarized Experiment
#' @examples 
#' gene.exp <- S4Vectors::DataFrame(sample1.exp = c("ENSG00000141510"=2.3,"ENSG00000171862"=5.4),
#'                   sample2.exp = c("ENSG00000141510"=1.6,"ENSG00000171862"=2.3))
#' dna.met <- S4Vectors::DataFrame(sample1.met = c("cg14324200"=0.5,"cg23867494"=0.1),
#'                        sample2.met =  c("cg14324200"=0.3,"cg23867494"=0.9))
#' createTSVTemplates(met = dna.met, exp = gene.exp)                       
#' @importFrom readr write_tsv
#' @export
createTSVTemplates <- function(met, exp) {
  assay <- c(rep("DNA methylation", ncol(met)),
             rep("Gene expression", ncol(exp)))
  primary <- rep("SampleX", ncol(met) + ncol(exp))
  colname <- c(colnames(met),colnames(exp))
  sampleMap <- data.frame(assay,primary,colname)
  message("===== Sample mapping example file ======")
  message("Saving example file as elmer_example_sample_mapping.tsv.")
  message("Please, fill primary column correctly")
  write_tsv(sampleMap,path = "elmer_example_sample_mapping.tsv")
  
  colData <- data.frame(primary = paste0("sample",1:ncol(met)), group = rep("To be filled",ncol(met)))
  message("===== Sample information example file ======")
  message("Saving example file as elmer_example_sample_metadata.tsv.")
  message("Please, fill primary column correctly, also you can add new columns as the example group column.")
  write_tsv(colData,path = "elmer_example_sample_metadata.tsv")
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


#' lable linear regression formula 
#' @param df A data.frame object contains two variables: dependent 
#' variable (Dep) and explanation variable (Exp).
#' @param Dep A character specify dependent variable. The first column 
#' will be dependent variable as default.
#' @param Exp A character specify explanation variable. The second column 
#' will be explanation variable as default.
#' @return A linear regression formula
#' @importFrom stats coef lm
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
#' # get GENCODE gene annotation (transcripts level)
#' \dontrun{
#'     getTSS <- getTSS()
#'     getTSS <- getTSS(genome.build = "hg38", TSS=list(upstream=1000, downstream=1000))
#' }
#' @export
#' @author Lijing Yao (maintainer: lijingya@usc.edu)
#' @import GenomeInfoDb
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom biomaRt useEnsembl
getTSS <- function(genome = "hg38",
                   TSS = list(upstream = NULL, downstream = NULL)){
  tries <- 0L
  msg <- character()
  while (tries < 3L) {
    tss <- tryCatch({
      if (genome == "hg19"){
        # for hg19
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           host = "feb2014.archive.ensembl.org",
                           path = "/biomart/martservice" ,
                           dataset = "hsapiens_gene_ensembl")
        attributes <- c("chromosome_name",
                        "start_position",
                        "end_position", 
                        "strand",
                        "transcript_start",
                        "transcript_end",
                        "ensembl_transcript_id",
                        "ensembl_gene_id", "entrezgene",
                        "external_gene_id")
      } else {
        # for hg38
        ensembl <- tryCatch({
          useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
        },  error = function(e) {
          useEnsembl("ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     mirror = "uswest")
        })
        attributes <- c("chromosome_name",
                        "start_position",
                        "end_position", "strand",
                        "ensembl_gene_id", 
                        "transcription_start_site",
                        "transcript_start",
                        "ensembl_transcript_id",
                        "transcript_end",
                        "external_gene_name")
      }
      chrom <- c(1:22, "X", "Y","M","*")
      db.datasets <- listDatasets(ensembl)
      description <- db.datasets[db.datasets$dataset=="hsapiens_gene_ensembl",]$description
      message(paste0("Downloading transcripts information. Using: ", description))
      
      filename <-  paste0(gsub("[[:punct:]]| ", "_",description),"_tss.rda")
      if(!file.exists(filename)) {
        tss <- getBM(attributes = attributes, filters = c("chromosome_name"), values = list(chrom), mart = ensembl)
        tss <- tss[!duplicated(tss$ensembl_transcript_id),]
        save(tss, file = filename)
      } else {
        tss <- get(load(filename))  
      } 
      if(genome == "hg19") tss$external_gene_name <- tss$external_gene_id
      tss$chromosome_name <-  paste0("chr", tss$chromosome_name)
      tss$strand[tss$strand == 1] <- "+" 
      tss$strand[tss$strand == -1] <- "-" 
      tss <- makeGRangesFromDataFrame(tss,start.field = "transcript_start", end.field = "transcript_end", keep.extra.columns = TRUE)
      
      if(!is.null(TSS$upstream) & !is.null(TSS$downstream)) 
        tss <- promoters(tss, upstream = TSS$upstream, downstream = TSS$downstream)
      tss
    }, error = function(e) {
      msg <<- conditionMessage(e)
      tries <<- tries + 1L
    })
    if(!is.null(tss)) break
  }
  if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)
  
  return(tss)
}

#' @title  Get human TF list from the UNiprot database
#' @description This function gets the last version of human TF list from the UNiprot database
#' @importFrom readr read_tsv
#' @return A data frame with the ensemble gene id.
getTF <- function(){
  print.header("Downloading TF list from Lambert, Samuel A., et al. The human transcription factors. Cell 172.4 (2018): 650-665.","subsection")
  # human.TF <- readr::read_table("http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt",col_names = F)
  # colnames(human.TF) <- "ensembl_gene_id"
  human.TF <- getdata("human.TF")
  return(human.TF)
}

#' @importFrom biomaRt getBM useMart listDatasets
get.GRCh <- function(genome = "hg19", genes, as.granges = FALSE) {
  tries <- 0L
  msg <- character()
  while (tries < 3L) {
    gene.location <- tryCatch({
      if (genome == "hg19"){
        # for hg19
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           host = "feb2014.archive.ensembl.org",
                           path = "/biomart/martservice" ,
                           dataset = "hsapiens_gene_ensembl")
        attributes <- c("ensembl_gene_id", "entrezgene","external_gene_id")
      } else {
        # for hg38
        ensembl <- tryCatch({
          useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
        },  error = function(e) {
          useEnsembl("ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     mirror = "uswest")
        })
        attributes <- c("ensembl_gene_id", "entrezgene","external_gene_name")
      }
      gene.location <- getBM(attributes = attributes,
                             filters = c("entrezgene"),
                             values = list(genes), mart = ensembl)
      colnames(gene.location) <-  c("ensembl_gene_id", "entrezgene","external_gene_name")
      gene.location <- gene.location[match(genes,gene.location$entrezgene),]
      gene.location
    }, error = function(e) {
      msg <<- conditionMessage(e)
      tries <<- tries + 1L
    })
    if(!is.null(gene.location)) break
  }
  if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)
  
  if(as.granges) {
    gene.location$strand[gene.location$strand == 1] <- "+"
    gene.location$strand[gene.location$strand == -1] <- "-"
    gene.location$chromosome_name <- paste0("chr",gene.location$chromosome_name)
    gene.location <- makeGRangesFromDataFrame(gene.location, seqnames.field = "chromosome_name",
                                              start.field = "start_position",
                                              end.field = "end_position",
                                              keep.extra.columns = TRUE) # considering the whole gene no their promoters
  }
  return(gene.location)
}



#' @title Get family of transcription factors
#' @description 
#' This will output a list each TF motif and TFs that binding the motis. Multiple TFs may
#' recognize a same motif such as TF family.  
#' The association between each motif famil and transcription factor was created using the 
#' (HOCOMOCO)[http://hocomoco.autosome.ru/human/mono] which TF structural families 
#' was created according to TFClass [@wingender2014tfclass]
#' This data is stored as a list whose elements 
#' are motifs and contents for each element are TFs which recognize the same motif that
#' is the name of the element. This data is used in function get.TFs in \pkg{ELMER} 
#' to identify the real regulator TF whose motif is enriched in a given set of probes 
#' and expression associate with average DNA methylation of these motif sites.
#' @importFrom rvest html_table  %>%
#' @importFrom xml2 read_html 
#' @export
#' @param classification Select if we will use Family classification or sub-family
#' @return A list of TFs and its family members
createMotifRelevantTfs <- function(classification = "family"){
  message("Retrieving TFClass ", classification," classification from ELMER.data.")
  if(classification == "family") motif.relevant.TFs <- getdata("TF.family")
  if(classification == "subfamily") motif.relevant.TFs <- getdata("TF.subfamily")
  return(motif.relevant.TFs)
}

#' @title Filtering probes
#' @description 
#' This function has some filters to the DNA methylation data
#' in each it selects probes to avoid correlations due to non-cancer 
#' contamination and for additional stringency.
#'  \itemize{
#' \item Filter 1: We usually call locus unmethylated when the methylation value < 0.3 and methylated when the methylation value > 0.3. 
#'       Therefore Meth_B is the percentage of methylation value > K. 
#'       Basically, this step will make sure we have at least a percentage of beta values lesser than K and n percentage of beta values greater K. 
#'       For example, if percentage is 5\%, the number of samples 100 and K = 0.3, 
#'       this filter will select probes that we have at least 5 (5\% of 100\%) samples have beta values > 0.3 and at least 5 samples have beta values < 0.3.
#'       This filter is importante as true promoters and enhancers usually have a pretty low value (of course purity can screw that up).  
#'       we often see lots of PMD probes across the genome with intermediate values like 0.4.  
#'       Choosing a value of 0.3 will certainly give some false negatives, but not compared to the number of false positives we thought we might get without this filter.
#' }
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' Method section (Linking enhancer probes with methylation changes to target genes with expression changes).
#' @param data A MultiAssayExperiment with a DNA methylation martrix or a DNA methylation matrix
#' @param K Cut off to consider probes as methylated or unmethylated. Default: 0.3
#' @param percentage The percentage of samples we should have at least considered as methylated and unmethylated
#' @return An object with the same class, but with the probes removed.
#' @importFrom MultiAssayExperiment experiments<-
#' @export
#' @examples
#'  random.probe <- runif(100, 0, 1)
#'  bias_l.probe <- runif(100, 0, 0.3)
#'  bias_g.probe <- runif(100, 0.3, 1)
#'  met <- rbind(random.probe,bias_l.probe,bias_g.probe)
#'  met <- preAssociationProbeFiltering(data = met,  K = 0.3, percentage = 0.05)
#'  met <- rbind(random.probe,random.probe,random.probe)
#'  met <- preAssociationProbeFiltering(met,  K = 0.3, percentage = 0.05)
#'  data <- ELMER:::getdata("elmer.data.example") # Get data from ELMER.data
#'  data <- preAssociationProbeFiltering(data,  K = 0.3, percentage = 0.05)
#'  
#'  cg24741609 <- runif(100, 0, 1)
#'  cg17468663 <- runif(100, 0, 0.3)
#'  cg14036402 <- runif(100, 0.3, 1)
#'  met <- rbind(cg24741609,cg14036402,cg17468663)
#'  colnames(met) <- paste("sample",1:100)
#'  exp <- met
#'  rownames(exp) <- c("ENSG00000141510","ENSG00000171862","ENSG00000171863")
#'  sample.info <- S4Vectors::DataFrame(sample.type = rep(c("Normal", "Tumor"),50))
#'  rownames(sample.info) <- colnames(exp)
#'  mae <- createMAE(exp = exp, met = met, colData = sample.info, genome = "hg38") 
#'  mae <- preAssociationProbeFiltering(mae,  K = 0.3, percentage = 0.05)
preAssociationProbeFiltering <- function(data, K = 0.3, percentage = 0.05){
  print.header("Filtering probes", type = "section")
  message("For more information see function preAssociationProbeFiltering")
  
  if(class(data) == class(MultiAssayExperiment())) { 
    met <- assay(getMet(data))
  } else {
    met <- data
  }
  # In percentage how many probes are bigger than K ?
  Meth_B <- rowMeans(met > K, na.rm = TRUE)
  # We should  have at least 5% methylation value < K or at least 5% methylation value > K
  keep <- Meth_B < (1 - percentage) & Meth_B > percentage
  message("Making sure we have at least ", percentage * 100, "% of beta values lesser than ", K," and ", 
          percentage * 100, "% of beta values greater ",K,".") 
  if(length(keep) - sum(keep) != 0) {
    message("Removing ", length(keep) - sum(keep), " probes out of ", length(keep))
  } else {
    message("There were no probes to be removed")
  }
  if(class(data) == class(MultiAssayExperiment())) {
    experiments(data)["DNA methylation"][[1]] <- experiments(data)["DNA methylation"][[1]][keep,]
  } else {
    data <- data[keep,,drop = FALSE]
  }
  return(data)
}

#' @importFrom xml2 read_html
#' @importFrom rvest html_table
getHocomocoTable <- function(){
  hocomoco <- tryCatch({
    hocomoco <- "http://hocomoco11.autosome.ru/human/mono?full=true" %>% read_html()  %>%  html_table()
    hocomoco <- hocomoco[[1]]
  }, error = function(e) {
    getdata("hocomoco.table")
  })
  
  TF.family <-  createMotifRelevantTfs()
  TF.subfamily <-  createMotifRelevantTfs("subfamily")
  
  x <- do.call(rbind.data.frame,lapply(TF.family, function(x) paste(x,collapse = ";")))
  x$Model <- names(TF.family)
  colnames(x) <- c("TF.family.member","Model")
  hocomoco <- merge(hocomoco,x, by = "Model")
  x <- do.call(rbind.data.frame,lapply(TF.subfamily, function(x) paste(x,collapse = ";")))
  x$Model <- names(TF.subfamily)
  colnames(x) <-  c("TF.subfamily.member","Model")
  hocomoco <- merge(hocomoco,x, by = "Model")
  return(hocomoco)
}

#' @title Get random pairs
#' @description 
#' This function will receive a pair gene probes and will return a 
#' random object with the following pattern, if a probe is linked to R1 and L3 genes
#' the random pairs will be a random probes (a distal probe not in the input pairs) 
#' also linked to its R1 and L3 gene.
#' @param pairs A data frame with probe, gene and side information. See example below.
#' @param met.platform DNA methyaltion platform to retrieve data from: EPIC or 450K (default)
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @param cores A interger which defines the number of cores to be used in parallel 
#' process. Default is 1: no parallel process.
#' @return A data frame with the random linkages
#' @export
#' @examples
#' \dontrun{
#'  data <- ELMER:::getdata("elmer.data.example")
#'  nearGenes <- GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
#'                             geneAnnot=getExp(data))
#'                             
#'  pair <- get.pair(data = data,
#'                   group.col = "definition", 
#'                   group1 = "Primary solid Tumor", 
#'                   group2 = "Solid Tissue Normal",
#'                   mode = "supervised",
#'                   diff.dir = "hypo",
#'                   nearGenes = nearGenes,
#'                   permu.size = 5,
#'                   raw.pvalue =  0.001,
#'                   Pe = 0.2,
#'                   dir.out="./",
#'                   permu.dir = "permu_test",
#'                   label = "hypo")
#' }
#'  pair <- data.frame(Probe = rep("cg00329272",3), 
#'                     GeneID = c("ENSG00000116213","ENSG00000130762","ENSG00000149527"),
#'                     Sides = c("R5","R2","L4"))                    
#'  getRandomPairs(pair)
getRandomPairs <- function(pairs, 
                           genome = "hg38",
                           met.platform = "450K",
                           cores = 1) {
  
  if(missing(pairs)) stop("Please set pairs argument")
  if(is.data.frame(pairs))
    pairs <- as.data.frame(pairs)
  
  # Get Probe information
  data("hm450.hg38.manifest")
  probes.ranges <- as.data.frame(hm450.hg38.manifest)[,1:5]
  colnames(probes.ranges) <- paste0("probe_",colnames(probes.ranges))
  
  # Get distal probes not in the pairs
  distal.probe <- get.feature.probe(genome = genome,
                                    met.platform = met.platform,
                                    feature = NULL) # get distal probes
  distal.probe <- distal.probe[!names(distal.probe) %in% pairs$Probe,] # Select probes were not used
  
  nb.pairs <- nrow(pairs)
  nb.probes <- length(unique(pairs$Probe))
  
  # get gene information
  genes <- TCGAbiolinks:::get.GRCh.bioMart(genome = genome,as.granges = TRUE)
  
  df.random <- NULL
  # We will get the double of random probes, because some will not be used in case it does not matches the position
  # Example: real probe + gene R10 and random probe does not have R10. Discart and get next random
  random.probes <- distal.probe[sample(1:length(distal.probe), nb.probes * 2),]
  near.genes <- GetNearGenes(TRange = random.probes, 
                             geneAnnot = genes, 
                             numFlankingGenes = 24, 
                             cores = cores)
  near.genes.df <- data.table::rbindlist(near.genes)
  # Now we should get the exactly same genes positions
  # if probe 1 was linked to R4 and L10 the random probe 1 will also be linked to its R4 and L10
  near.genes.linked <- NULL
  eval <- 1
  for(p in 1:length(near.genes)){
    side <- unique(pairs[pairs$Probe == unique(pairs$Probe)[eval],"Sides"])
    aux <-  near.genes.df[near.genes.df$Target == names(near.genes)[p],]
    same <- aux[aux$Side %in% side,]
    # If I do not have the same nearby positions
    if(length(side) != nrow(same)) next
    near.genes.linked <- rbind(near.genes.linked, same)
    eval <- eval+ 1
    if(length(unique(near.genes.linked$Target)) == nb.probes) break
  }
  colnames(near.genes.linked)[1] <- "Probe" 
  
  # Add probe metadata to output
  probes.ranges$Probe <- rownames(probes.ranges)
  near.genes.linked <- merge(near.genes.linked,probes.ranges,by ="Probe")
  
  return(near.genes.linked)
}


# Reading homer output. For each reagion (rows)
# homer will try to find if a given motif was found in it.
# This will read this homer file and create a sparce matrix
# in which 1 means the motif was found and 0 not found.
# this is used to calculate the motif enrichement compared 
# to a background signal.
getMatrix <- function(filename) {
  motifs <- readr::read_tsv(filename)
  # From 1 to 21 we have annotations
  matrix <- Matrix::Matrix(0, nrow = nrow(motifs), ncol = ncol(motifs) - 21 ,sparse = TRUE)
  colnames(matrix) <- gsub(" Distance From Peak\\(sequence,strand,conservation\\)","",colnames(motifs)[-c(1:21)])
  rownames(matrix) <- motifs$PeakID
  matrix[!is.na(motifs[,-c(1:21)])] <- 1
  matrix <- as(matrix, "nsparseMatrix")
  return(matrix)
}

#' @title Calculate motif Erichment 
#' @description Calculates fisher exact test
#' @param foreground A nsparseMatrix object in each 1 means the motif is found in a region, 0 not.
#' @param background A nsparseMatrix object in each 1 means the motif is found in a region, 0 not.
#' @examples 
#' foreground <- Matrix::Matrix(sample(0:1,size = 100,replace = T), nrow = 10, ncol = 10,sparse = TRUE)
#' rownames(foreground) <- paste0("region",1:10)
#' colnames(foreground) <- paste0("motif",1:10)
#' background <- Matrix::Matrix(sample(0:1,size = 100,replace = T), nrow = 10, ncol = 10,sparse = TRUE)
#' rownames(background) <- paste0("region",1:10)
#' colnames(background) <- paste0("motif",1:10)
#' calculateEnrichement(foreground,background)
calculateEnrichement <- function(foreground, 
                                 background){
  if(missing(foreground)) stop("foreground argument is missing")
  if(missing(background)) stop("background argument is missing")
  
  # a is the number of probes within the selected probe set that contain one or more motif occurrences; 
  # b is the number of probes within the selected probe set that do not contain a motif occurrence; 
  # c and d are the same counts within the entire enhancer probe set (background)
  # lower boundary of 95% conf idence interval = exp (ln OR - SD)
  a <- Matrix::colSums(foreground)
  b <- nrow(foreground) - Matrix::colSums(foreground)
  c <- Matrix::colSums(background)
  d <- nrow(background) - Matrix::colSums(background)
  fisher <- plyr::adply(seq_len(length(a)),.margins = 1, .fun = function(i)  { 
    x <- fisher.test(matrix(c(a[i],b[i],c[i],d[i]),nrow = 2,ncol = 2))
    ret <- data.frame(x$conf.int[1],x$conf.int[2],x$estimate,x$p.value)
    colnames(ret) <- c("lowerOR","upperOR","OR","p.value")
    ret
  },.id = NULL,.progress = "text")
  rownames(fisher) <- names(a)
  Summary <- data.frame(motif  =  names(a),
                        NumOfRegions = Matrix::colSums(foreground, na.rm=TRUE),
                        fisher,
                        FDR = p.adjust(fisher$p.value,method = "BH"),
                        stringsAsFactors = FALSE)
  Summary <- Summary[order(-Summary$lowerOR),]
  return(Summary)
}


#' @title Use Hocomoco motif and homer to identify motifs in a given region
#' @param regions A GRanges object. Names will be used as the identifier.
#' @param output.filename Final file name
#' @param region.size If NULL the motif will be mapped to the region. If set a window around its center will be considered.
#' For example if region.size is 500, then +-250bp round it will be searched.
#' @param cores A interger which defines the number of cores to be used in parallel 
#' process. Default is 1: no parallel process.
#' @param genome Homer genome (hg38, hg19)
#' @description 
#' To find for each probe the know motif we will use HOMER software (http://homer.salk.edu/homer/).
#' Homer and genome should be installed before this function is executed
#' Step:
#' 1 - get DNA methylation probes annotation with the regions
#' 2 - Make a bed file from it
#' 3 - Execute section: Finding Instance of Specific Motifs from http://homer.salk.edu/homer/ngs/peakMotifs.html to the HOCOMOCO TF motifs
#' Also, As HOMER is using more RAM than the available we will split the files in to 100k probes.
#' Obs: for each probe we create a winddow of 500 bp (-size 500) around it. This might lead to false positives, but will not have false negatives.
#' The false posives will be removed latter with some statistical tests.
#' @examples 
#' \dontrun{
#'  # use the center of the region and +-250bp around it
#'  gr0 <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), IRanges(1:10, width=10:1))
#'  names(gr0) <- paste0("ID",c(1:10))
#'  findMotifRegion(regions = gr0, region.size = 500, genome = "hg38", cores = 1)
#'  
#'  # use the region size itself
#'  gr1 <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), IRanges(1:10, width=sample(200:1000,10)))
#'  names(gr1) <- paste0("ID",c(1:10))
#'  findMotifRegion(regions = gr0, genome = "hg38", cores = 1)
#' }
findMotifRegion <- function(regions, 
                            output.filename = "mapped_motifs_regions.txt",
                            region.size = NULL,
                            genome = "hg38",
                            cores = 1){
  
  if(!is(regions, class(GRanges()))) stop("Regions must be a Genomic Ranges object")
  # get all hocomoco 11 motifs
  TFBS.motif <- "http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif"
  if(!file.exists(basename(TFBS.motif))) downloader::download(TFBS.motif,basename(TFBS.motif))
  
  if(!file.exists(output.filename)){
    df <- data.frame(seqnames = seqnames(gr),
                     starts = as.integer(start(gr)),
                     ends = end(gr),
                     names = names(gr),
                     scores = c(rep(".", length(gr))),
                     strands = strand(gr))
    n <- nrow(df)
    step <- ifelse(n > 1000, 1000, n )
    pb <- txtProgressBar(max = floor(n/step), style = 3)
    
    for(j in 0:floor(n/step)){
      setTxtProgressBar(pb, j)
      # STEP 2
      file.aux <- paste0(gsub(".txt","",output.filename),"_",j,".bed")
      if(!file.exists(gsub(".bed",".txt",file.aux))){
        end <- ifelse(((j + 1) * step) > n, n,((j + 1) * step))
        write.table(df[((j * step) + 1):end,], file = file.aux, col.names = F, quote = F,row.names = F,sep = "\t")
        
        # STEP 3 use -mscore to get scores
        # we need to check if annotatePeaks.pl is available!
        cmd <- paste("annotatePeaks.pl" ,file.aux, genome, "-m", basename(TFBS.motif), 
                     ifelse(is.null(region.size),"",paste0("-size ", region.size)),
                     "-cpu", cores,
                     ">", gsub(".bed",".txt",file.aux))
        error <- system(cmd)
        if(error == 127) {
          unlink(file.aux)
          unlink(gsub(".bed",".txt",file.aux))
          stop("annotatePeaks.pl had an error. Please check homer install")
        }
      }
    }
  }
  close(pb)
  # We will merge the results from each file into one
  peaks <- NULL
  pb <- txtProgressBar(max = floor(n/step), style = 3)
  for(j in 0:floor(n/step)){
    setTxtProgressBar(pb, j)
    f <- paste0(gsub(".txt","",output.filename),"_",j,".txt")
    if(file.exists(f)){
      aux <-  readr::read_tsv(f)
      colnames(aux)[1] <- "PeakID"
      if(is.null(peaks)) {
        peaks <- aux
      } else {
        peaks <- rbind(peaks, aux)
      }
    }
    gc()
  }
  close(pb)
  print("Writing file...")
  if(!is.null(peaks)) readr::write_tsv(peaks, path = file,col_names = TRUE)
  print("DONE!")
}
