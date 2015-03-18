## Anaylsis TCGA data pipeline------------
#' ELMER analysis pipe for TCGA data.
#' @param disease TCGA short form disease name such as COAD
#' @param analysis a vector of characters listing the analysis need to be done. 
#' Analysis are "download","distal.enhancer","diffMeth","pair","motif","TF.search". 
#' Default is "all" meaning all the analysis will be processed. 
#' @param wd a path showing working dirctory. Default is "./"
#' @param cores A interger which defines number of core to be used in parallel process. 
#' Default is NULL: don't use parallel process.
#' @param Data A path showing the folder containing DNA methylation, expression and clinic data
#' @param ... A list of parameters for functions: GetNearGenes, get.feature.probe, 
#' get.diff.meth, get.pair,
#' @return Different analysis results.
#' @export 
#' @examples
#' distal.probe <- TCGA.pipe(disease = "LUSC", analysis="distal.enhancer", wd="~/")
TCGA.pipe <- function(disease,analysis="all",wd="./",cores=NULL,Data=NULL,...){
  if(missing(disease)) 
    stop("Disease should be specified.\nDisease short name (such as LAML) 
         please check https://tcga-data.nci.nih.gov/tcga/.")
  if(analysis[1] == "all") analysis=c("download","distal.enhancer",
                                      "diffMeth","pair","motif","TF.search")
  disease <- toupper(disease)
  dir.out <- sprintf("%s/Result/%s",wd,disease)
  if(!file.exists(dir.out)) dir.create(dir.out,recursive = TRUE)
  args <- list(...)
  #Download 
  if("download" %in% analysis){
    message("###################\nDownload data\n###################\n\n")
    if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    params <- args[names(args) %in% c("RNAtype","Methfilter")]
    params$disease <- disease
    params$basedir <- sprintf("%s/Data",wd)
    do.call(getTCGA,params)
    analysis <- analysis[!analysis %in% "download"]
  }
    
  ## select distal enhancer probes
  if("distal.enhancer" %in% analysis){
    message("###################\nSelect distal enhancer probes\n###################\n\n")
    params <- args[names(args) %in% c("probe","TSS","feature","TSS.range","distal")]
    probeInfo <- do.call(get.feature.probe,params)
    save(probeInfo,file=sprintf("%s/probeInfo_feature.rda",dir.out))
    if(length(analysis)==1){
      return(probeInfo)
    }else{
      analysis <- analysis[!analysis %in% "distal.enhancer"]
      invisible(gc())
    }
  }
  
  #get differential DNA methylation
  if("diffMeth" %in% analysis){
    message("###################\nGet differential DNA methylation loci\n###################\n\n")
    #refine meth data: filter out non tumor or normal samples--------------
    meth.file <- sprintf("%s/%s_meth_refined.rda",dir.out,disease)
    if(!file.exists(meth.file)){
      if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
      load(sprintf("%s/probeInfo_feature.rda",dir.out))
      load(sprintf("%s/%s_meth.rda",Data,disease))
      TN <- sapply(colnames(Meth),tcgaSampleType)
      Meth <- Meth[rownames(Meth) %in% as.character(probeInfo$name),TN %in% c("Tumor","Normal")]
      save(Meth,file= sprintf("%s/%s_meth_refined.rda",dir.out,disease))
      rm(Meth)
    }
    mee <- fetch.mee(meth=meth.file,TCGA=TRUE,
                     probeInfo=sprintf("%s/probeInfo_feature.rda",dir.out))
    params <- args[names(args) %in% c("diff.dir","percentage","pvalue","sig.dif")]
    params <- c(params,list(dir.out=dir.out, cores=cores))
    diff.meth <- do.call(get.diff.meth,c(params,list(mee=mee)))
    if(length(analysis)==1){
      return(diff.meth)
    }else{
      analysis <- analysis[!analysis %in% "diffMeth"]
      rm(mee,diff.meth)
      invisible(gc())
    }
  }
  
  
  #predict pair
  if("pair" %in% analysis){
    message("###################\nPredict pairs\n###################\n\n")
    meth.file <- sprintf("%s/%s_meth_refined.rda",dir.out,disease)
    exp.file <- sprintf("%s/%s_RNA_refined.rda",dir.out,disease)
    if(!file.exists(exp.file)){
      load(sprintf("%s/%s_RNA.rda",Data,disease))
      TN <- sapply(colnames(GeneExp),tcgaSampleType)
      GeneExp <- GeneExp[,TN %in% c("Tumor","Normal")]
      GeneExp <- log2(GeneExp+1)
      save(GeneExp,file= sprintf("%s/%s_RNA_refined.rda",dir.out,disease))
      rm(GeneExp)
    }
    ## construct geneAnnot for finding nearby gene
    geneAnnot <- args[names(args) %in% "geneAnnot"]
    if(length(geneAnnot)==0){
      geneAnnot <- sprintf("%s/geneInfo.rda",dir.out)
      if(!file.exists(geneAnnot)){
        newenv <- new.env()
        load(system.file("extdata","UCSC_gene_hg19.rda",package = "ELMER"),
             envir=newenv)
        txs <- get(ls(newenv)[1],envir=newenv)
        geneInfo <- promoters(txs,upstream = 0, downstream = 0)
        geneInfo$GENEID <- paste0("ID",geneInfo$GENEID)
        save(geneInfo,file=geneAnnot)
      }
    }else{
      geneAnnot <- geneAnnot[["geneAnnot"]]
    }
    distal.probe <- suppressWarnings(get.feature.probe(feature="distal"))
    mee<- fetch.mee(meth=meth.file, exp=exp.file, probeInfo=distal.probe, 
                    geneInfo=geneAnnot,TCGA=TRUE)
    ## define diff.dir
    diff.dir <- args[names(args) %in% "diff.dir"]
    if(length(diff.dir)==0 || diff.dir[["diff.dir"]]=="both"){
      diff.dir <- c("hyper","hypo")
    }else{
      diff.dir <- diff.dir[["diff.dir"]]
    }
    ## calculation
    SigPair <- list()
    for(diff in diff.dir){
      message(sprintf("Identify putative probe-gene pair for %smethylated probes",diff))
      #Construct data.
      Sig.probes <- read.csv(sprintf("%s/getMethdiff.%s.probes.significant.csv",dir.out,diff),
                             stringsAsFactors=FALSE)[,1]
      ## Get nearby genes-----------------------
      nearGenes.file <- args[names(args) %in% "nearGenes"]
      if(length(nearGenes.file)==0){
        nearGenes.file <- sprintf("%s/%s.probes_nearGenes.rda",dir.out,diff)
        if(!file.exists(nearGenes.file)){
          params <- args[names(args) %in% c("geneNum")]
          nearGenes <- do.call(GetNearGenes,
                               c(list(TRange=getProbeInfo(mee,probe=Sig.probes),
                                      geneAnnot=getGeneInfo(mee),cores=cores),
                                 params))
          save(nearGenes,file=nearGenes.file)
        }
      }else{
        nearGenes.file <- nearGenes.file[["nearGenes"]]
      }
      ## get pair
      permu.dir <- paste0(dir.out,"/permu")
      params <- args[names(args) %in% c("percentage","permu.size","Pe")]
      SigPair[[diff]] <- do.call(get.pair,c(list(mee=mee,probes=Sig.probes,
                                                 nearGenes=nearGenes.file,
                                                 permu.dir=permu.dir,
                                                 dir.out=dir.out,
                                                 cores=cores,label=diff),
                                            params))
    }
    if(length(analysis)==1){
      return(SigPair)
    }else{
      analysis <- analysis[!analysis %in% "pair"]
      rm(mee,SigPair)
      invisible(gc())
    }
  }
  
  # search enriched motif
  if("motif" %in% analysis){
    message("###################\nMotif search\n###################\n\n")
    ## define diff.dir
    diff.dir <- args[names(args) %in% "diff.dir"]
    if(length(diff.dir)==0 || diff.dir[["diff.dir"]]=="both"){
      diff.dir <- c("hyper","hypo")
    }else{
      diff.dir <- diff.dir[["diff.dir"]]
    }
    for(diff in diff.dir){
      message(sprintf("Identify enriched motif for %smethylated probes",diff))
      if(file.exists(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, diff))){
        Sig.probes <- unique(read.csv(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, diff), 
                                      stringsAsFactors=FALSE)$Probe)
      }else{
        stop(sprintf("%s/%s.pairs.significant.csv file doesn't exist",dir.out, diff))
      }
      params <- args[names(args) %in% c("background.probes","lower.OR","min.incidence")]
      enriched.motif <- do.call(get.enriched.motif, c(list(probes=Sig.probes,
                                                           dir.out=dir.out,
                                                           label=diff),
                                                      params))
    }
    if(length(analysis)==1){
      return(enriched.motif)
    }else{
      analysis <- analysis[!analysis %in% "motif"]
      rm(enriched.motif)
      invisible(gc())
    }
  }
   
   #search responsible TFs
  if("TF.search" %in% analysis){
    message("###################\nSearch responsible TFs\n###################\n\n")
    ## make mee 
    #construct RNA seq data
    meth.file <- sprintf("%s/%s_meth_refined.rda",dir.out,disease)
    exp.file <- sprintf("%s/%s_RNA_refined.rda",dir.out,disease)
    probeInfo <- sprintf("%s/probeInfo_feature.rda",dir.out)
    geneAnnot <- sprintf("%s/geneInfo.rda",dir.out)
    mee <- fetch.mee(meth=meth.file,exp=exp.file,TCGA=TRUE,
                     probeInfo=probeInfo,geneInfo=geneAnnot)
    diff.dir <- args[names(args) %in% "diff.dir"]
    if(length(diff.dir)==0 || diff.dir[["diff.dir"]]=="both"){
      diff.dir <- c("hyper","hypo")
    }else{
      diff.dir <- diff.dir[["diff.dir"]]
    }
    for(diff in diff.dir){
      message(sprintf("Identify regulatory TF for enriched motif in %smethylated probes",diff))
      enriched.motif <- args[names(args) %in% "enriched.motif"]
      if(length(enriched.motif)==0){
        enriched.motif <- sprintf("%s/getMotif.%s.enriched.motifs.rda",dir.out, diff)
      }
      params <- args[names(args) %in% c("TFs", "motif.relavent.TFs","percentage")]
      do.call(get.TFs, c(list(mee=mee, enriched.motif=enriched.motif,
                              dir.out=dir.out, cores=cores, 
                              label=diff),
                         params))
    }
  }
}

