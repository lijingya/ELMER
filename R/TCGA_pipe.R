## Anaylsis TCGA data pipeline------------
#' ELMER analysis pipe for TCGA data.
#' @param disease TCGA short form disease name such as COAD
#' @param analysis a vector of characters listing the analysis need to be done. Analysis are "download","distal.enhancer","diffMeth","pair","motif","TF.search". Default is "all" meaning all the analysis will be processed. 
#' @param wd a path showing working dirctory. Default is "./"
#' @param cores A interger which defines number of core to be used in parallel process. Default is NULL: don't use parallel process.
#' @param Data A path showing the folder containing DNA methylation, expression and clinic data
#' @param ... A list of parameters for functions: GetNearGenes, get.feature.probe, get.diff.meth, get.pair,
#' @return Different analysis results.
#' @export 

TCGA.pipe <- function(disease,analysis="all",wd="./",cores=NULL,Data=NULL,...){
  if(missing(disease)) stop("Disease should be specified.\nDisease short name (such as LAML) please check https://tcga-data.nci.nih.gov/tcga/.")
  if(analysis == "all") analysis=c("download","distal.enhancer","diffMeth","pair","motif","TF.search")
  disease <- toupper(disease)
  if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
  dir.out <- sprintf("%s/Result/%s",wd,disease)
  if(!file.exists(dir.out)) dir.create(dir.out,recursive = T)
  args <- list(...)
  browser()
  #Download 
  if("download" %in% analysis){
    cat("########################################\nDownload data\n########################################")
    params <- args[names(args) %in% c("RNAtype","Methfilter")]
    params$disease <- disease
    params$basedir <- sprintf("%s/Data",wd)
    do.call(getTCGA,params)
    analysis <- analysis[!analysis %in% "download"]
  }else if(any(grepl("rda",dir(Data)))){
    stop("Data need to be specified.\nEither provide downloaded data folder contain your data stored as XXX.rda.")
  }
    
  ## select distal enhancer probes
  if("distal.enhancer" %in% analysis){
    cat("########################################\nSelect distal enhancer probes\n########################################")
    params <- args[names(args) %in% c("probe","TSS","feature","TSS.range","distal")]
    probeInfo <- do.call(get.feature.probe,params)
    save(probeInfo,file=sprintf("%s/probeInfo_feature.rda",dir.out))
    if(length(analysis)==1){
      return(probeInfo)
    }else{
      analysis <- analysis[!analysis %in% "distal.enhancer"]
    }
  }
  
  #get differential DNA methylation
  if("diffMeth" %in% analysis){
    cat("########################################\nGet differential DNA methylation loci\n########################################")
    #refine meth data: filter out non tumor or normal samples--------------
    meth.file <- sprintf("%s/%s_meth_refined.rda",dir.out,disease)
    if(!file.exists(meth.file)){
      load(sprintf("%s/probeInfo_feature.rda",dir.out))
      load(sprintf("%s/%s_meth.rda",Data,disease))
      TN <- sapply(colnames(Meth),tcgaSampleType)
      Meth <- Meth[rownames(Meth) %in% as.character(probeInfo$name),TN %in% c("Tumor","Normal")]
      save(Meth,file= sprintf("%s/%s_meth_refined.rda",dir.out,disease))
      rm(Meth)
    }
    mee <- fetch.data(meth=meth.file,TCGA=T,probeInfo=sprintf("%s/probeInfo_feature.rda",dir.out))
    params <- args[names(args) %in% c("diff.dir","percentage","pvalue","sig.dif")]
    params <- c(params,list(dir.out=dir.out, cores=cores))
    diff.meth <- do.call(get.diff.meth,c(params,list(mee=mee)))
    if(length(analysis)==1){
      return(diff.meth)
    }else{
      analysis <- analysis[!analysis %in% "diffMeth"]
    }
  }
  
  #predict pair
  if("pair" %in% analysis){
    cat("########################################\nPredict pairs\n########################################")
    #construct RNA seq data
    meth.file <- sprintf("%s/%s_meth_refined.rda",dir.out,disease)
    exp.file <- sprintf("%s/%s_RNA_refined.rda",dir.out,disease)
    probeInfo <- sprintf("%s/probeInfo_feature.rda",dir.out)
    if(!file.exists(exp.file)){
      load(sprintf("%s/%s_RNA.rda",Data,disease))
      TN <- sapply(colnames(GeneExp),tcgaSampleType)
      GeneExp <- GeneExp[,TN %in% c("Tumor","Normal")]
      save(GeneExp,file= sprintf("%s/%s_RNA_refined.rda",dir.out,disease))
      rm(GeneExp)
    }
    ## construct geneAnnot for finding nearby gene
    geneAnnot <- args[names(args) %in% "geneAnnot"]
    if(length(geneAnnot)==0){
      geneAnnot <- sprintf("%s/geneInfo.rda",dir.out)
      if(!file.exists(geneAnnot)){
        newenv <- new.env()
        load(system.file("extdata","UCSC_gene_hg19.rda",package = "ELMER"),envir=newenv)
        txs <- get(ls(newenv)[1],envir=newenv)
        geneInfo <- promoters(txs,upstream = 0, downstream = 0)
        save(geneInfo,file=geneAnnot)
      }
    }else{
      geneAnnot <- geneAnnot[["geneAnnot"]]
    }
    
    ## define diff.dir
    diff.dir <- args[names(args) %in% "diff.dir"]
    if(length(diff.dir)==0 || diff.dir[["diff.dir"]]=="both"){
      diff.dir <- c("hyper","hypo")
    }else{
      diff.dir <- diff.dir[["diff.dir"]]
    }
    
    ## calculation
    for(diff in diff.dir){
      #Construct data.
      Sig.probes <- read.csv(sprintf("%s/%s.probes_significant.csv",dir.out,diff),stringsAsFactors=F)[,1]
      mee <- fetch.data(meth=meth.file,exp=exp.file,TCGA=T,probeInfo=probeInfo,geneInfo=geneAnnot)
      ## Get nearby genes-----------------------
      nearGenes.file <- args[names(args) %in% "nearGenes"]
      if(length(nearGenes.file)==0){
        nearGenes.file <- sprintf("%s/%s.probes_nearGenes.rda",dir.out,diff)
        if(!file.exists(nearGenes.file)){
          params <- args[names(args) %in% c("geneNum")]
          nearGenes <- do.call(GetNearGenes,c(list(TRange=getProbeInfo(mee,probe=Sig.probes),geneAnnot=mee@geneInfo,cores=cores),params))
          save(nearGenes,file=nearGenes.file)
        }
      }else{
        nearGenes.file <- nearGenes.file[["nearGenes"]]
      }
      ## get pair
      permu.dir <- paste0(dir.out,"/permu")
      #(mee,probes,nearGenes,percentage=0.2,permu.size=1000,permu.dir=NULL, Pe=0.01,dir.out="./",cores=4,label=NULL)
      params <- args[names(args) %in% c("percentage","permu.size","Pe")]
      do.call(get.pair,c(list(mee=mee,probes=Sig.probes,nearGenes=nearGenes.file,permu.dir=permu.dir,dir.out=dir.out,cores=cores,label=diff)))
    }
  }
  
  # search enriched motif
  if("motif" %in% analysis){
    cat("########################################\nMotif search\n########################################")
    ## define diff.dir
    diff.dir <- args[names(args) %in% "diff.dir"]
    if(length(diff.dir)==0 || diff.dir[["diff.dir"]]=="both"){
      diff.dir <- c("hyper","hypo")
    }else{
      diff.dir <- diff.dir[["diff.dir"]]
    }
    for(diff in diff.dir){
      if(file.exists(sprintf("%s/%s.significant.pairs.csv",dir.out, diff))){
        Sig.probes <- read.csv(sprintf("%s/%s.significant.pairs.csv",dir.out, label), stringsAsFactors=F)$Probe
      }else{
        stop(sprintf("%s/%s.significant.pairs.csv file doesn't exist",dir.out, diff))
      }
      params <- args[names(args) %in% c("background.probes","lower.OR","min.incidence")]
      enriched.motif <- do.call(get.enriched.motif, c(list(probes=Sig.probes,dir.out=dir.out,label=diff),params))
    }
   
   
  }
#   
#   #search responsible TFs
#   if("TF.search" %in% analysis){
#     cat("########################################\nSearch responsible TFs\n########################################")
# #    params <- args[names(args) %in% c(need to be done)]
# #    TF.search
#   }
}

