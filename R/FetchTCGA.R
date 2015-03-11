#' getTCGA
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param Meth A logic if TRUE HM450K DNA methylation data will download.
#' @param RNA A logic if TRUE RNA-seq Hiseq-V2 from TCGA level 3 will be download.
#' @param Clinic A logic if TRUE clinic data will be download for that disease.
#' @param basedir A path shows where the data will be stored.
#' @param RNAtype A charactor to specify whether use isoform level or gene level. When RNAtype=gene, gene level gene expression will be used. When isoform, then isoform data will be used.
#' @param Methfilter A number. For each probe, the percentage of NA among the all the samples should smaller than Methfilter.
#' @export
getTCGA <- function(disease,Meth=TRUE,RNA=TRUE,Clinic=TRUE,basedir="./Data",RNAtype="gene",Methfilter=0.2){
  if(missing(disease)) stop("disease need to be specified.")
  if(Meth){
    message("########################################\nDownloading methylation\n########################################\n\n")
    test.meth <- tryCatch({get450K(disease, basedir)},error=function(err){return("error")})
    if(!test.meth == "error") matrixMeth(disease,basedir,filter=Methfilter)
  }
  
  if(RNA){
    message("########################################\nDownloading RNA\n########################################\n\n")
    test.rna <- tryCatch({getRNAseq(disease, basedir)},error=function(err){return("error")})
    if(!test.rna == "error") matrixRNA(disease,basedir,RNAtype)
  }
 
  if(Clinic){
    message("########################################\nDownloading Clinic \n########################################\n\n")
    test.clinic <- tryCatch({getClinic(disease, basedir)},error=function(err){return("error")})
    if(!test.clinic == "error") matrixClinic(disease,basedir)
  }
  if(test.meth=="error") warning(sprintf("Failed to download DNA methylation data. Possible possibility: 1. No 450K DNA methylation data for %s patients; 2.Wget doesn't work.",disease))
  if(test.rna=="error") warning(sprintf("Failed to download RNA-seq data. Possible possibility: 1. No RNA-seq data for %s patients; 2.Wget doesn't work.",disease))
  if(test.clinic=="error") warning(sprintf("Failed to download clinic data. Possible possibility:1. No clinical data for %s patients; 2.Wget doesn't work.",disease))
}

#' getRNAseq
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @export
getRNAseq <- function(disease, basedir="./Data")
{
  wd <- getwd()
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.rna <- file.path(dir.raw,"RNA")
  if(!file.exists(dir.rna)) dir.create(dir.rna,recursive = T)
  setwd(dir.rna)
  target <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor"
  domain <- "cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2"
  baseURL <- paste(target, disease, domain, sep="/")
  links <- system(sprintf("wget --no-check-certificate -O - %s | grep -o '<a href=['\"'\"'\"][^\"'\"'\"']*['\"'\"'\"]' | sed -e 's/^<a href=[\"'\"'\"']//' -e 's/[\"'\"'\"']$//' | grep -v 'tar'",baseURL), intern=TRUE)
  mage <- links[grep("mage-tab", links)]
  RNA <- links[grep("RNASeqV2.Level_3", links)]
  mage.versions <- unlist(lapply(strsplit(mage,"[.]"),function(x){return(as.numeric(x[[6]]))}))
  RNA.versions <- unlist(lapply(strsplit(RNA,"[.]"),function(x){return(as.numeric(x[[6]]))}))
  mage.version <- sub("/", ".tar.gz",mage[mage.versions == max(mage.versions)])
  RNA.version <- sub("/", ".tar.gz",RNA[RNA.versions == max(RNA.versions)])
  mage.URL <- paste(baseURL, mage.version, sep="/")
  RNA.URL <- paste(baseURL, RNA.version, sep="/")
  if(length(mage.URL)==0 | length(RNA.URL)==0){
    setwd(wd)
    stop("No RNA-seq data.")
  }else{
    cmd <- list(mage = paste("wget --no-check-certificate -O",mage.version, mage.URL, "-a", "mage_download.log", sep=" "), 
                RNA = paste("wget --no-check-certificate -O",RNA.version, RNA.URL, "-a",  "RNA_download.log", sep=" "))
    lapply(cmd, system)
    system(sprintf("tar -zxvf %s", mage.version),ignore.stdout=T)
    system(sprintf("tar -zxvf %s", RNA.version),ignore.stdout=T)
    system (sprintf("rm %s %s", mage.version, RNA.version))
    setwd(wd)
  } 
}

#' get450K
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @export
get450K <- function(disease, basedir="./Data")
{
  wd <- getwd()
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.meth <- file.path(dir.raw,"Meth")
  if(!file.exists(dir.meth)) dir.create(dir.meth,recursive = T)
  setwd(dir.meth)
  target <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor"
  domain <- "cgcc/jhu-usc.edu/humanmethylation450/methylation/"
  baseURL <- paste(target, disease, domain, sep="/")
  links <- system(sprintf("wget --no-check-certificate -O - %s | grep -o '<a href=['\"'\"'\"][^\"'\"'\"']*['\"'\"'\"]' | sed -e 's/^<a href=[\"'\"'\"']//' -e 's/[\"'\"'\"']$//' | grep -v 'tar'",baseURL), intern=TRUE)
  mage <- links[grep("mage-tab", links)]
  Meth <- links[grep("HumanMethylation450.Level_3", links)]
  mage.versions <- unlist(lapply(strsplit(mage,"[.]"),function(x){return(as.numeric(x[[6]]))}))
  meth.versions <- unlist(lapply(strsplit(Meth,"[.]"),function(x){return(as.numeric(x[[6]]))}))
  mage.version <- sub("/", ".tar.gz",mage[mage.versions == max(mage.versions)])
  meth.version <- sub("/", ".tar.gz",Meth[meth.versions == max(meth.versions)])
  mage.URL <- paste(baseURL, mage.version, sep="/")
  Meth.URL <- paste(baseURL,meth.version, sep="/")
  
  if(length(mage.URL)==0 | length(mage.URL)==0){
    setwd(wd)
    stop("No 450K DNA methylation data.")
  }else{
    cmd <- list(mage = paste("wget --no-check-certificate -O",mage.version, mage.URL, "-a", "mage_download.log", sep=" "), 
                Meth = paste(paste("wget --no-check-certificate -O",meth.version, Meth.URL, "-a",  "Meth_download.log", sep=" "),collapse = ";"))
    lapply(cmd, system)
    system(sprintf("tar -zxvf %s", mage.version),ignore.stdout=T)
    system(paste(sprintf("tar -zxvf %s", meth.version),collapse=";"),ignore.stdout=T)
    system (sprintf("rm %s %s", mage.version, paste(meth.version,collapse = " ")))
    setwd(wd)
  }
}

#' getClinic
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @export
getClinic <- function(disease, basedir="./Data")
{
  wd <- getwd()
  disease <- tolower(disease)
  diseasedir <- file.path(basedir, toupper(disease))
  dir.raw <- file.path(diseasedir,"Raw")
  dir.clinic <- file.path(dir.raw,"Clinic")
  if(!file.exists(dir.clinic)) dir.create(dir.clinic,recursive = T)
  setwd(dir.clinic)
  target <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor"
  domain <- "bcr/biotab/clin/"
  baseURL <- paste(target, disease, domain, sep="/")
  links <- system(sprintf("wget --no-check-certificate -O - %s | grep -o '<a href=['\"'\"'\"][^\"'\"'\"']*['\"'\"'\"]' | sed -e 's/^<a href=[\"'\"'\"']//' -e 's/[\"'\"'\"']$//' | grep -v 'tar'",baseURL), intern=TRUE)
  clinic <- links[grep("clinical_patient", links)]
  if(length(clinic)==0){
    setwd(wd)
    stop(sprintf("No clinical data for %s patients.",disease))
  }else{
    Clinic.URL <- paste(baseURL,clinic, sep="/")
    cmd <-paste("wget --no-check-certificate -O",clinic, Clinic.URL, "-a", "download.log", sep=" ")
    system(cmd)
    setwd(wd)
  }
}

## type is gene level expression, isoform level expression
#' matrixRNA
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @param type A charactor to specify whether use isoform level or gene level. When RNAtype=gene, gene level gene expression will be used. When isoform, then isoform data will be used.
#' @export
matrixRNA <- function(disease,basedir="./Data",type="gene"){
  if(missing(disease)) stop("disease should be specified.")
  wd <- getwd()
  disease <- tolower(disease)
  diseasedir <- file.path(wd,basedir, toupper(disease))
  dir.RNA <- file.path(diseasedir,"Raw","RNA")
  ## read file names
  mage.file <- dir(dir(path = dir.RNA,pattern = "mage-tab",full.names = T),"sdrf.txt", full.names = T)
  Files <- read.table(mage.file, sep="\t", stringsAsFactors = F, header=T)[,c(2,22)]
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
  GeneExp <- do.call(cbind,mclapply(Files, function(x){read.table(x,stringsAsFactors = F, sep="\t",header=T)[,2]},mc.cores=detectCores()/2))
  GENEID <- read.table(Files[1], stringsAsFactors = F, sep="\t", header=T)[,1]
  rownames(GeneExp) <- GENEID
  GeneExp <- GeneIDName(GeneExp)
  save(GeneExp,file=sprintf("%s/%s_RNA.rda",diseasedir,toupper(disease)))
}

#' matrixMeth
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @param filter For each probe, the percentage of NA among the all the samples should smaller than filter.
#' @export
matrixMeth <- function(disease,basedir="./Data",filter=0.2){
  if(missing(disease)) stop("disease should be specified.")
  wd <- getwd()
  disease <- tolower(disease)
  diseasedir <- file.path(wd,basedir, toupper(disease))
  dir.meth <- file.path(diseasedir,"Raw","Meth")
  ## read file names
  mage.file <- dir(dir(path = dir.meth,pattern = "mage-tab",full.names = T),"sdrf.txt", full.names = T)
  Files <- unique(read.table(mage.file, sep="\t", stringsAsFactors = F, header=T)[,c(2,29,28)])
  ID <- Files[,1]
  Files <- paste0(dir.meth,"/",paste(Files[,2],Files[,3],sep="/"))
  names(Files) <- ID
  ##load data
  Meth <- do.call(cbind,mclapply(Files, function(x){read.table(x,stringsAsFactors = F, sep="\t",header=T,skip=1)[,2]},mc.cores=detectCores()/2))
  Probe <- read.table(Files[1], stringsAsFactors = F, sep="\t", header=T,skip=1)[,1]
  rownames(Meth) <- Probe
  Meth <- Meth[rowMeans(is.na(Meth))<0.2,]
  save(Meth,file=sprintf("%s/%s_meth.rda",diseasedir,toupper(disease)))
}


#' matrixClinic
#' @param disease A character to specify disease in TCGA such as BLCA
#' @param basedir A path shows where the data will be stored.
#' @export
matrixClinic <- function(disease,basedir="./Data"){
  if(missing(disease)) stop("disease should be specified.")
  wd <- getwd()
  disease <- tolower(disease)
  diseasedir <- file.path(wd,basedir, toupper(disease))
  dir.Clinic <- file.path(diseasedir,"Raw","Clinic")
  File <- dir(dir.Clinic,pattern = "nationwidechildrens.org_clinical_patient",full.names = T)
  Clinic <- read.table(File,sep="\t",header=T,stringsAsFactors = F)[-c(1:2),]
  Useful <- unlist(apply(Clinic,2,function(x){tmp <- unique(x) 
                                          if(length(tmp)<2 & all(grepl("Not",tmp))){return(FALSE)
                                                                                    }else{
                                                                                      return(TRUE)
                                                                                    }}))
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
