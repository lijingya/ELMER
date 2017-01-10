#' getTCGA to download DNA methylation, RNA expression and clinic data for all samples of certain cancer type from TCGA.
#' @description 
#' getTCGA is a function to download DNA methylation, RNA expression and clinic data for all
#'  samples of certain cancer type from TCGA website. And downloaded data will be transform
#'   to matrixes or data frame for further analysis.
#' @param disease A character specifies the disease to download in TCGA such as BLCA
#' @param Meth A logic if TRUE HM450K DNA methylation data will download.
#' @param RNA A logic if TRUE RNA-seq Hiseq-V2 from TCGA level 3 will be download.
#' @param Clinic A logic if TRUE clinic data will be download for that disease.
#' @param genome Data aligned against which genome of reference. Options: "hg19", "hg38" (default)
#' @param basedir A path shows where the data will be stored.
#' @param RNAtype A charactor to specify whether use isoform level or gene level. 
#' When RNAtype=gene, gene level gene expression will be used.
#' @param Methfilter A number. For each probe, the percentage of NA among the all the samples should smaller than Methfilter.
#' @return Download DNA methylation (HM450K)/RNAseq(HiseqV2)/Clinic data for
#' the specified disease from TCGA.
#' @usage getTCGA(disease, Meth = TRUE, RNA = TRUE, Clinic = TRUE, 
#'                basedir = "./Data", RNAtype = "gene", Methfilter = 0.2)
#' @export
#' @examples
#' getTCGA("BRCA",Meth=FALSE, RNA=FALSE, Clinic=TRUE, basedir="~")
getTCGA <- function(disease,
                    Meth=TRUE,
                    RNA=TRUE,
                    Clinic=TRUE,
                    basedir="./Data",
                    RNAtype="gene",
                    genome = "hg38",
                    Methfilter=0.2){
  if(missing(disease)) stop("disease need to be specified.")
  if(Meth){
    print.header("Downloading DNA methylation", "subsection")
    test.meth <- tryCatch({
      get450K(disease, basedir,filter = Methfilter, genome = genome)
    }, error = function(err){
      return("error")
    })
  }
  
  if(RNA){
    print.header("Downloading RNA", "subsection")
    test.rna <- tryCatch({
      getRNAseq(disease, basedir, genome = genome)
    }, error=function(err){
      return("error")
    })
  }
  
  if(Clinic){
    print.header("Downloading Clinic", "subsection")
    test.clinic <- tryCatch({
      getClinic(disease, basedir)
    }, error=function(err){
      return("error")
    })
  }
  
  if(Meth && test.meth == "error") 
    warning(
      sprintf("Failed to download DNA methylation data. Possible possibility: 
              1. No 450K DNA methylation data for %s patients; 
              2. Download error.",disease))
  if(RNA && test.rna == "error") 
    warning(
      sprintf("Failed to download RNA-seq data. Possible possibility: 
               1. No RNA-seq data for %s patients; 
               2. Download error.",disease))
  if(Clinic && test.clinic == "error") 
    warning(
      sprintf("Failed to download clinic data. Possible possibility:
               1. No clinical data for %s patients; 
               2. Download error.", disease))
}

#' getRNAseq to download all RNAseq data for a certain cancer type from TCGA.
#' @description getRNAseq is a function to download RNAseq data for all samples of a certain cancer type from TCGA
#' @param disease A character specifies disease in TCGA such as BLCA
#' @param basedir Download all RNA seq level 3 data for the specified disease.
#' @param genome Data aligned against which genome of reference. Options: "hg19", "hg38" (default)
#' @usage getRNAseq(disease, basedir = "./Data")
#' @return Download all RNA seq level 3 data for the specified disease.
#' @importFrom TCGAbiolinks GDCdownload GDCquery GDCprepare
getRNAseq <- function(disease, 
                      basedir="./Data",
                      genome = "hg38") {
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.rna <- file.path(dir.raw,"RNA")
  if (!file.exists(dir.rna)) dir.create(dir.rna,recursive = TRUE)
  
  if (genome == "hg38"){
    query <- GDCquery(project = paste0("TCGA-",toupper(disease)), 
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM-UQ")
  } else {
    query <- GDCquery(project = paste0("TCGA-",toupper(disease)),
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type  = "normalized_results",
                      experimental.strategy = "RNA-Seq",
                      legacy = TRUE)
  }
  tryCatch({
    GDCdownload(query, directory = dir.rna, chunks.per.download = 200)
  }, error = function(e) {
    GDCdownload(query, directory = dir.rna, chunks.per.download = 50)
  })
  
  # Preparing to save output if it does not exists
  fout <- sprintf("%s/%s_RNA.rda",diseasedir,toupper(disease))
  if(!file.exists(fout)){
    rna <- GDCprepare(query, 
                      directory = dir.rna,
                      summarizedExperiment = TRUE)
    if(genome == "hg19"){
      rownames(rna) <- values(rna)$ensembl_gene_id
    }
    message(paste0("Saving Gene Expression to: ", fout))
    save(rna,file=fout)
  } else {
    message(paste("Gene Expression object already exists:", fout))
  }
  return("OK")
}

#' get450K to download HM40K DNA methylation data for certain cancer types from TCGA website.
#'  @description 
#'  get450K is a function to download latest version of HM450K DNA methylation 
#'  for  all samples of certain cancer types from GDC website.
#' @param disease A character specifies the disease to download from TCGA such as BLCA
#' @param basedir A path. Shows where the data will be stored.
#' @param filter For each probe, the percentage of NA among the all the samples
#'  should smaller than filter.
#' @param genome Data aligned against which genome of reference. Options: "hg19", "hg38" (default)
#' @return Download all DNA methylation from HM450K level 3 data for
#'  the specified disease.
#' @importFrom TCGAbiolinks GDCquery GDCdownload GDCprepare
#' @usage get450K(disease, basedir = "./Data")
get450K <- function(disease, 
                    basedir="./Data",
                    filter=0.2,
                    genome = "hg38"){
  
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.meth <- file.path(dir.raw,"Meth")
  if (!file.exists(dir.meth)) dir.create(dir.meth,recursive = TRUE)
  
  if (genome == "hg38") {
    query <- GDCquery(project = paste0("TCGA-",toupper(disease)), 
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450")
  } else {
    query <- GDCquery(project = paste0("TCGA-",toupper(disease)), 
                      data.category = "DNA methylation",
                      legacy = TRUE,
                      platform = "Illumina Human Methylation 450")
  }  
  tryCatch({
    GDCdownload(query,directory = dir.meth,chunks.per.download = 5)
  }, error = function(e) {
    GDCdownload(query,directory = dir.meth,  method = "client")
  })
  
  # Preparing to save output if it does not exists
  fout <- sprintf("%s/%s_meth.rda",diseasedir,toupper(disease))
  if(!file.exists(fout)){
    met <- GDCprepare(query = query,
                      directory = dir.meth,
                      summarizedExperiment = TRUE)
    
    # Remove probes that has more than 20% of its values as NA
    met <- met[rowMeans(is.na(assay(met))) < filter,]
    message(paste0("Saving DNA methylation to: ", fout))
    save(met,file = fout)
  } else {
    message(paste("DNA methylation object already exists:", fout))
  }
  return("OK")
}

#' getClinic to download clinic data for certain cancer types from TCGA website.
#' @description 
#' getClinic is a function to download latest version of clinic data for all samples of certain cancer types from TCGA website.
#' @param disease A character specifies the disease to download from TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @importFrom TCGAbiolinks GDCquery_clinic
#' @return Download all clinic information for the specified disease.
getClinic <- function(disease, basedir="./Data")
{
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.clinic <- file.path(dir.raw,"Clinic")
  if(!file.exists(dir.clinic)) dir.create(dir.clinic,recursive = TRUE)
  
  Clinic <- GDCquery_clinic(project = paste0("TCGA-",toupper(disease)))
  save(Clinic,file=sprintf("%s/%s_clinic.rda",diseasedir,toupper(disease)))
  return("OK")
}

print.header <- function(text, type ="section"){
  message(paste(rep("-",nchar(text) + 3),collapse = ""))
  message(paste(ifelse(type=="section","o","oo"),text))
  message(paste(rep("-",nchar(text)+ 3),collapse = ""))
}
