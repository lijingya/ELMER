#' getTCGA to download DNA methylation, RNA expression and clinic data for all samples of certain cancer type from TCGA.
#' @description 
#' getTCGA is a function to download DNA methylation, RNA expression and clinic data for all
#'  samples of certain cancer type from TCGA website. And downloaded data will be transform
#'   to matrixes or data frame for further analysis.
#' @param disease A character specifies the disease to download in TCGA such as BLCA
#' @param Meth A logic if TRUE HM450K DNA methylation data will download.
#' @param RNA A logic if TRUE RNA-seq Hiseq-V2 from TCGA level 3 will be download.
#' @param Clinic A logic if TRUE clinic data will be download for that disease.
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
                    Methfilter=0.2){
  if(missing(disease)) stop("disease need to be specified.")
  if(Meth){
    message("################\nDownloading DNA methylation\n################\n\n")
    test.meth <- tryCatch({get450K(disease, basedir,filter = Methfilter)},
                          error=function(err){return("error")})
  }
  
  if(RNA){
    message("################\nDownloading RNA\n################\n\n")
    test.rna <- tryCatch({getRNAseq(disease, basedir)},
                         error=function(err){return("error")})
  }
  
  if(Clinic){
    message("################\nDownloading Clinic \n################\n\n")
    test.clinic <- tryCatch({getClinic(disease, basedir)},
                            error=function(err){return("error")})
  }
  
  if(Meth && test.meth=="error") 
    warning(
      sprintf("Failed to download DNA methylation data. Possible possibility: 
              1. No 450K DNA methylation data for %s patients; 
              2.Wget doesn't work.",disease))
  if(RNA && test.rna=="error") 
    warning(
      sprintf("Failed to download RNA-seq data. Possible possibility: 
              1. No RNA-seq data for %s patients; 
              2.Wget doesn't work.",disease))
  if(Clinic && test.clinic=="error") 
    warning(
      sprintf("Failed to download clinic data. Possible possibility:
              1. No clinical data for %s patients; 
              2.Wget doesn't work.",disease))
}

#' getRNAseq to download all RNAseq data for a certain cancer type from TCGA.
#' @description getRNAseq is a function to download RNAseq data for all samples of a certain cancer type from TCGA
#' @param disease A character specifies disease in TCGA such as BLCA
#' @param basedir Download all RNA seq level 3 data for the specified disease.
#' @usage getRNAseq(disease, basedir = "./Data")
#' @return Download all RNA seq level 3 data for the specified disease.
getRNAseq <- function(disease, basedir="./Data")
{
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.rna <- file.path(dir.raw,"RNA")
  if(!file.exists(dir.rna)) dir.create(dir.rna,recursive = TRUE)
  query <- GDCquery(project = paste0("TCGA-",toupper(disease)), 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", file.type  = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  GDCdownload(query)
  rna <- GDCprepare(query)
  GeneExp <- TCGAprepare_elmer(rna,platform = "IlluminaHiSeq_RNASeqV2", save = TRUE)
  save(GeneExp,file=sprintf("%s/%s_RNA.rda",diseasedir,toupper(disease)))
}

#' get450K to download HM40K DNA methylation data for certain cancer types from TCGA website.
#'  @description 
#'  get450K is a function to download latest version of HM450K DNA methylation 
#'  for  all samples of certain cancer types from GDC website.
#' @param disease A character specifies the disease to download from TCGA such as BLCA
#' @param basedir A path. Shows where the data will be stored.
#' @param filter For each probe, the percentage of NA among the all the samples
#'  should smaller than filter.
#' @return Download all DNA methylation from HM450K level 3 data for
#'  the specified disease.
#' @importFrom TCGAbiolinks GDCquery GDCdownload GDCprepare
#' @usage get450K(disease, basedir = "./Data")
get450K <- function(disease, 
                    basedir="./Data",
                    filter=0.2){
  
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.meth <- file.path(dir.raw,"Meth")
  if(!file.exists(dir.meth)) dir.create(dir.meth,recursive = TRUE)
  
  query <- GDCquery(project = paste0("TCGA-",toupper(disease)), 
                    legacy = TRUE,
                    data.category = "DNA methylation",
                    platform = "Illumina Human Methylation 450")
  
  json  <- tryCatch(GDCdownload(query,directory = dir.meth, chunks.per.download = 10),
                    error = function(e) GDCdownload(query, method = "client"))
  met <- GDCprepare(query = query,
                    directory = dir.meth,
                    save = TRUE, 
                    save.filename = paste(disease,"_DNA_methylation_450K.rda"),
                    summarizedExperiment = TRUE)
  Meth <- TCGAprepare_elmer(rna,platform = "IlluminaHiSeq_RNASeqV2", met.na.cut = 0.2, save = TRUE)
  save(Meth,file=sprintf("%s/%s_meth.rda",diseasedir,toupper(disease)))
}

#' getClinic to download clinic data for certain cancer types from TCGA website.
#' @description 
#' getClinic is a function to download latest version of clinic data for all samples of certain cancer types from TCGA website.
#' @param disease A character specifies the disease to download from TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @return Download all clinic information for the specified disease.
getClinic <- function(disease, basedir="./Data")
{
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.clinic <- file.path(dir.raw,"Clinic")
  if(!file.exists(dir.clinic)) dir.create(dir.clinic,recursive = TRUE)
  
  clin.query <- GDCquery(project = paste0("TCGA-",toupper(disease)), 
                         data.category = "Clinical", barcode = "TCGA-FD-A5C0")
  json  <- tryCatch(GDCdownload(clin.query, directory = dir.clinic), 
                    error = function(e) GDCdownload(clin.query, method = "client",directory = dir.clinic))
  clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient",directory = dir.clinic)
  save(Clinic,file=sprintf("%s/%s_clinic.rda",diseasedir,toupper(disease)))
}

#Gene make rowname separat -------------------------------------
GeneIDName <- function(x){
  tmp<-strsplit(rownames(x),"\\|")
  GeneID<-unlist(lapply(tmp,function(x) x[2]))
  GeneID <- paste0("ID",GeneID)
  row.names(x) <- GeneID
  return(x)
}
