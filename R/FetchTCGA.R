#' getTCGA
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param Meth A logic if TRUE HM450K DNA methylation data will download.
#' @param RNA A logic if TRUE RNA-seq Hiseq-V2 from TCGA level 3 will be download.
#' @param Clinic A logic if TRUE clinic data will be download for that disease.
#' @param basedir A path shows where the data will be stored.
#' @param RNAtype A charactor to specify whether use isoform level or gene level.
#'  When RNAtype=gene, gene level gene expression will be used. When isoform, 
#'  then isoform data will be used.
#' @param Methfilter A number. For each probe, the percentage of NA among the 
#' all the samples should smaller than Methfilter.
#' @return Download DNA methylation (HM450K)/RNAseq(HiseqV2)/Clinic data for
#' the specified disease from TCGA.
#' @export
#' @examples
#' getTCGA("BRCA",Meth=FALSE, RNA=FALSE, Clinic=TRUE, basedir="~")
getTCGA <- function(disease,Meth=TRUE,RNA=TRUE,Clinic=TRUE,basedir="./Data",
                    RNAtype="gene",Methfilter=0.2){
  if(missing(disease)) stop("disease need to be specified.")
  if(Meth){
    message("################\nDownloading methylation\n################\n\n")
    test.meth <- tryCatch({get450K(disease, basedir)},
                          error=function(err){return("error")})
    if(!test.meth == "error") matrixMeth(disease,basedir,filter=Methfilter)
  }
  
  if(RNA){
    message("################\nDownloading RNA\n################\n\n")
    test.rna <- tryCatch({getRNAseq(disease, basedir)},
                         error=function(err){return("error")})
    if(!test.rna == "error") matrixRNA(disease,basedir,RNAtype)
  }
 
  if(Clinic){
    message("################\nDownloading Clinic \n################\n\n")
    test.clinic <- tryCatch({getClinic(disease, basedir)},
                            error=function(err){return("error")})
    if(!test.clinic == "error") matrixClinic(disease,basedir)
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

#' getRNAseq
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @return Download all RNA seq level 3 data for the specified disease.
getRNAseq <- function(disease, basedir="./Data")
{
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.rna <- file.path(dir.raw,"RNA")
  if(!file.exists(dir.rna)) dir.create(dir.rna,recursive = TRUE)
  target <- 
    "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor"
  domain <- "cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2"
  baseURL <- paste(target, disease, domain, sep="/")
  cmd <- paste0("wget --no-check-certificate -O - ",baseURL,
                " | grep -o '<a href=['\"'\"'\"][^\"'\"'\"']*['\"'\"'\"]' | ",
                "sed -e 's/^<a href=[\"'\"'\"']//' -e 's/[\"'\"'\"']$//'",
                " | grep -v 'tar'")
  links <-system(cmd, intern=TRUE)
  mage <- links[grep("mage-tab", links)]
  RNA <- links[grep("RNASeqV2.Level_3", links)]
  mage.versions <- unlist(lapply(strsplit(mage,"[.]"),
                                 function(x){return(as.numeric(x[[6]]))}))
  RNA.versions <- unlist(lapply(strsplit(RNA,"[.]"),
                                function(x){return(as.numeric(x[[6]]))}))
  mage.version <- sub("/", ".tar.gz",mage[mage.versions == max(mage.versions)])
  RNA.version <- sub("/", ".tar.gz",RNA[RNA.versions == max(RNA.versions)])
  mage.URL <- paste(baseURL, mage.version, sep="/")
  RNA.URL <- paste(baseURL, RNA.version, sep="/")
  if(length(mage.URL)==0 | length(RNA.URL)==0){
    stop("No RNA-seq data.")
  }else{
    cmd <- list(mage = paste("wget --no-check-certificate -O",
                             paste0(dir.rna,"/",mage.version), mage.URL, "-a", 
                             "mage_download.log", sep=" "), 
                RNA = paste("wget --no-check-certificate -O",
                            paste0(dir.rna,"/",RNA.version), RNA.URL, "-a",  
                            "RNA_download.log",sep=" "))
    lapply(cmd, system)
    system(sprintf("tar -zxvf %s", paste0(dir.rna,"/",mage.version)),ignore.stdout=TRUE)
    system(sprintf("tar -zxvf %s", paste0(dir.rna,"/",RNA.version)),ignore.stdout=TRUE)
    system (sprintf("rm %s %s", paste0(dir.rna,"/",mage.version), paste0(dir.rna,"/",RNA.version)))
  } 
}

#' get450K
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @return Download all DNA methylation from HM450K level 3 data for the
#' specified disease.
get450K <- function(disease, basedir="./Data")
{
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.meth <- file.path(dir.raw,"Meth")
  if(!file.exists(dir.meth)) dir.create(dir.meth,recursive = TRUE)
  target <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor"
  domain <- "cgcc/jhu-usc.edu/humanmethylation450/methylation/"
  baseURL <- paste(target, disease, domain, sep="/")
  cmd <- paste0("wget --no-check-certificate -O - ",baseURL,
                " | grep -o '<a href=['\"'\"'\"][^\"'\"'\"']*['\"'\"'\"]' ",
                "| sed -e 's/^<a href=[\"'\"'\"']//' -e 's/[\"'\"'\"']$//' ",
                "| grep -v 'tar'")
  links <- system(cmd, intern=TRUE)
  mage <- links[grep("mage-tab", links)]
  Meth <- links[grep("HumanMethylation450.Level_3", links)]
  mage.versions <- unlist(lapply(strsplit(mage,"[.]"),
                                 function(x){return(as.numeric(x[[6]]))}))
  meth.versions <- unlist(lapply(strsplit(Meth,"[.]"),
                                 function(x){return(as.numeric(x[[6]]))}))
  mage.version <- sub("/", ".tar.gz",mage[mage.versions == max(mage.versions)])
  meth.version <- sub("/", ".tar.gz",Meth[meth.versions == max(meth.versions)])
  mage.URL <- paste(baseURL, mage.version, sep="/")
  Meth.URL <- paste(baseURL,meth.version, sep="/")
  
  if(length(mage.URL)==0 | length(mage.URL)==0){
    stop("No 450K DNA methylation data.")
  }else{
    cmd <- list(mage = paste("wget --no-check-certificate -O",
                             paste0(dir.meth,"/",mage.version), mage.URL,"-a", 
                             "mage_download.log", sep=" "), 
                Meth = paste(paste("wget --no-check-certificate -O", 
                                   paste0(dir.meth,"/",meth.version), Meth.URL, "-a", 
                                   "Meth_download.log", sep=" "),
                             collapse = ";"))
    lapply(cmd, system)
    system(sprintf("tar -zxvf %s", paste0(dir.meth,"/",mage.version)),ignore.stdout=TRUE)
    system(paste(sprintf("tar -zxvf %s", paste0(dir.meth,"/",meth.version)),collapse=";"),
           ignore.stdout=TRUE)
    system (sprintf("rm %s %s", paste0(dir.meth,"/",mage.version), 
                    paste(paste0(dir.meth,"/",meth.version),collapse = " ")))
  }
}

#' getClinic
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @return Download all clinic information for the specified disease.
getClinic <- function(disease, basedir="./Data")
{
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.clinic <- file.path(dir.raw,"Clinic")
  if(!file.exists(dir.clinic)) dir.create(dir.clinic,recursive = TRUE)
  target <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor"
  domain <- "bcr/biotab/clin/"
  baseURL <- paste(target, disease, domain, sep="/")
  cmd <- paste0("wget --no-check-certificate -O - ",baseURL,
                " | grep -o '<a href=['\"'\"'\"][^\"'\"'\"']*['\"'\"'\"]' ",
                "| sed -e 's/^<a href=[\"'\"'\"']//' -e 's/[\"'\"'\"']$//' ",
                "| grep -v 'tar'")
  links <- system(cmd, intern=TRUE)
  clinic <- links[grep("clinical_patient", links)]
  if(length(clinic)==0){
    stop(sprintf("No clinical data for %s patients.",disease))
  }else{
    Clinic.URL <- paste(baseURL,clinic, sep="/")
    cmd <-paste("wget --no-check-certificate","-O",paste0(dir.clinic,"/",clinic), Clinic.URL, "-a", 
                "download.log", sep=" ")
    system(cmd)
  }
}

## type is gene level expression, isoform level expression
#' matrixRNA
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @param type A charactor to specify whether use isoform level or gene level.
#'  When RNAtype=gene, gene level gene expression will be used. When isoform, 
#'  then isoform data will be used.
#' @return A matrix of gene expression values for samples (column) and genes (row)
matrixRNA <- function(disease,basedir="./Data",type="gene"){
  if(missing(disease)) stop("disease should be specified.")
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.RNA <- file.path(diseasedir,"Raw","RNA")
  ## read file names
  mage.file <- dir(dir(path = dir.RNA,pattern = "mage-tab",full.names = TRUE),
                   "sdrf.txt", full.names = TRUE)
  Files <- read.table(mage.file, sep="\t", stringsAsFactors = FALSE, header=TRUE)[,c(2,22)]
  if(type %in% "isoform"){
    Files <- Files[grepl("rsem.isoforms.normalized_results",Files[,2]),]
  }else{
    Files <- Files[grepl("rsem.genes.normalized_results",Files[,2]),]
  }
  ID <- Files[,1]
  subfolder <- list.dirs(dir.RNA)[grepl("HiSeq_RNASeqV2.Level",list.dirs(dir.RNA))]
  Files <- paste0(subfolder,"/",Files[,2])
  names(Files) <- ID
  ##load data
  if(requireNamespace("parallel", quietly=TRUE)) {
  GeneExp <- do.call(cbind,parallel::mclapply(Files, 
                                    function(x){
                                      read.table(x,stringsAsFactors = FALSE, 
                                                 sep="\t",header=TRUE)[,2]},
                                    mc.cores=parallel::detectCores()/2))
  } else {
  GeneExp <- do.call(cbind,lapply(Files, 
								  function(x){
									  read.table(x,stringsAsFactors = FALSE, 
												 sep="\t",header=TRUE)[,2]}))
  }

  GENEID <- read.table(Files[1], stringsAsFactors = FALSE, sep="\t", header=TRUE)[,1]
  rownames(GeneExp) <- GENEID
  GeneExp <- GeneIDName(GeneExp)
  save(GeneExp,file=sprintf("%s/%s_RNA.rda",diseasedir,toupper(disease)))
}

#' matrixMeth
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @param filter For each probe, the percentage of NA among the all the samples
#'  should smaller than filter.
#' @return A matrix of DNA methylation values for samples (column) and probes (row)
matrixMeth <- function(disease,basedir="./Data",filter=0.2){
  if(missing(disease)) stop("disease should be specified.")
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.meth <- file.path(diseasedir,"Raw","Meth")
  ## read file names
  mage.file <- dir(dir(path = dir.meth,pattern = "mage-tab",full.names = TRUE),
                   "sdrf.txt", full.names = TRUE)
  Files <- unique(read.table(mage.file, sep="\t", stringsAsFactors = FALSE, header=TRUE)[,c(2,29,28)])
  ID <- Files[,1]
  Files <- paste0(dir.meth,"/",paste(Files[,2],Files[,3],sep="/"))
  names(Files) <- ID
  ##load data
  if(requireNamespace("parallel", quietly=TRUE)) {
	  Meth <- do.call(cbind,parallel::mclapply(Files, 
											   function(x){
												   read.table(x,stringsAsFactors = FALSE, sep="\t",header=TRUE,skip=1)[,2]},
												   mc.cores=parallel::detectCores()/2))
  } else {
	  Meth <- do.call(cbind,lapply(Files, 
								   function(x){
									   read.table(x,stringsAsFactors = FALSE, sep="\t",header=TRUE,skip=1)[,2]}))
  }
  Probe <- read.table(Files[1], stringsAsFactors = FALSE, sep="\t", header=TRUE,skip=1)[,1]
  rownames(Meth) <- Probe
  Meth <- Meth[rowMeans(is.na(Meth))<0.2,]
  save(Meth,file=sprintf("%s/%s_meth.rda",diseasedir,toupper(disease)))
}


#' matrixClinic
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @return A data frame contains clinic information for samples
matrixClinic <- function(disease,basedir="./Data"){
  if(missing(disease)) stop("disease should be specified.")
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.Clinic <- file.path(diseasedir,"Raw","Clinic")
  File <- dir(dir.Clinic,pattern = "nationwidechildrens.org_clinical_patient",full.names = TRUE)
  Clinic <- read.delim(File,sep="\t",header=TRUE,stringsAsFactors = FALSE)[-c(1:2),]
  Useful <- unlist(apply(Clinic,2,function(x){tmp <- unique(x) 
                                          if(length(tmp)<2 & all(grepl("Not",tmp)))
                                            {return(FALSE)
                                             }else{
                                             return(TRUE)}}))
  Clinic <- Clinic[,Useful]
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
