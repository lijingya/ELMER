#' ELMER analysis pipeline for TCGA data.
#' @description 
#' ELMER analysis pipeline for TCGA data. This pipeline combine every steps of \pkg{ELMER}
#' analyses: get.feature.probe, get.diff.meth, get.pair, get.permu, get.enriched.motif and get.TFs.
#' Every steps' results are saved.
#' @param disease TCGA short form disease name such as COAD
#' @param genome Data aligned against which genome of reference. Options: "hg19", "hg38" (default)
#' @param analysis A vector of characters listing the analysis need to be done.
#' Analysis can be "download","distal.probes","diffMeth","pair","motif","TF.search".
#' Default is "all" meaning all the analysis will be processed.
#' @param wd A path shows working dirctory. Default is "./"
#' @param cores A interger which defines number of core to be used in parallel process. 
#' Default is 1: don't use parallel process.
#' @param diff.dir A character can be "hypo" or "hyper", showing dirction DNA methylation changes.  
#' If it is "hypo", get.diff.meth function will identify all significantly hypomethylated
#' CpG sites; If "hyper", get.diff.meth function will identify all significantly hypermethylated
#' CpG sites
#' @param Data A path shows the folder containing DNA methylation, expression and clinic data
#' @param ... A list of parameters for functions: GetNearGenes, get.feature.probe, 
#' get.diff.meth, get.pair
#' @return Different analysis results.
#' @export 
#' @examples
#' \dontrun{
#'   distal.probe <- TCGA.pipe(disease = "LUSC", analysis="Probe.selection", wd="~/")
#' }
TCGA.pipe <- function(disease,
                      genome = "hg38",
                      analysis="all",
                      wd="./",
                      cores=1,
                      Data=NULL, 
                      diff.dir="hypo",
                      ...){
  if (missing(disease)) 
    stop("Disease should be specified.\nDisease short name (such as LAML) 
         please check https://gdc-portal.nci.nih.gov")
  
  available.analysis <- c("download","distal.enhancer",
                          "diffMeth","pair",
                          "motif","TF.search","all")
  if (analysis[1] == "all") analysis = available.analysis[1:6]
  
  if(any(!analysis %in% tolower(analysis))) 
    stop(paste0("Availbale options for analysis argument are: ",
                paste(c("",available.analysis), collapse = "\n=> ")))
  
  disease <- toupper(disease)
  dir.out <- sprintf("%s/Result/%s",wd,disease)
  if(!file.exists(dir.out)) dir.create(dir.out,recursive = TRUE)
  args <- list(...)
  
  # Download 
  if("download" %in% tolower(analysis)){
    message("###################\nDownload data\n###################\n\n")
    if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    params <- args[names(args) %in% c("RNAtype","Methfilter")]
    params$disease <- disease
    params$basedir <- sprintf("%s/Data",wd)
    do.call(getTCGA,params)
    analysis <- analysis[!analysis %in% "download"]
  }
  
  ## select distal enhancer probes
  if("distal.probes" %in% tolower(analysis)){
    message("###################\nSelect distal enhancer probes\n###################\n\n")
    params <- args[names(args) %in% c("TSS","feature","TSS.range","rm.chr")]
    probeInfo <- do.call(get.feature.probe,params)
    save(probeInfo,file = sprintf("%s/probeInfo_feature_distal.rda",dir.out))
    if(length(analysis) == 1){
      return(probeInfo)
    } else {
      analysis <- analysis[!analysis %in% "distal.probes"]
      invisible(gc())
    }
  }
  #get differential DNA methylation
  if("diffMeth" %in% tolower(analysis)){
    message("###################\nGet differential DNA methylation loci\n###################\n\n")
    #refine meth data: filter out non tumor or normal samples--------------
    meth.file <- sprintf("%s/%s_meth_refined.rda",dir.out,disease)
    if(!file.exists(meth.file)){
      if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
      load(sprintf("%s/%s_meth.rda",Data,disease))
      TN <- sapply(colnames(Meth),tcgaSampleType)
      Meth <- Meth[,TN %in% c("Tumor","Normal")]
      save(Meth,file = sprintf("%s/%s_meth_refined.rda",dir.out,disease))
      rm(Meth)
    }
    mae <- createMAE(met = meth.file, TCGA=TRUE, filter.probes = sprintf("%s/probeInfo_feature_distal.rda",dir.out))
    #mee <- fetch.mee(meth=meth.file,TCGA=TRUE,
    #                 probeInfo=sprintf("%s/probeInfo_feature_distal.rda",dir.out))
    
    params <- args[names(args) %in% c("percentage","pvalue","sig.dif")]
    params <- c(params,list(diff.dir=diff.dir, dir.out=dir.out, cores=cores))
    diff.meth <- do.call(get.diff.meth,c(params,list(data=mae)))
    if(length(analysis)==1){
      return(diff.meth)
    }else{
      analysis <- analysis[!analysis %in% "diffMeth"]
      rm(mae,diff.meth)
      invisible(gc())
    }
  }
  
  
  #predict pair
  if("pair" %in% tolower(analysis)){
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
    if(length(geneAnnot) == 0){
      geneAnnot <- sprintf("%s/geneInfo.rda",dir.out)
      if (!file.exists(geneAnnot)){
        geneInfo <- txs(TSS = list(upstream=0, downstream=0))
        geneInfo$GENEID <- paste0("ID",geneInfo$GENEID)
        save(geneInfo,file = geneAnnot)
      }
    } else {
      geneAnnot <- geneAnnot[["geneAnnot"]]
    }
    ## get distal probe info
    distal.probe <- sprintf("%s/probeInfo_feature_distal.rda",dir.out)
    if(!file.exists(distal.probe)){
      params <- args[names(args) %in% c("TSS","feature","TSS.range","rm.chr")]
      distal.probe <- suppressWarnings(do.call(get.feature.probe,params))
    }
    
    mee <- fetch.mee(meth=meth.file, exp=exp.file, probeInfo=distal.probe, 
                     geneInfo=geneAnnot,TCGA=TRUE)
    
    ## calculation
    message(sprintf("Identify putative probe-gene pair for %smethylated probes",diff.dir))
    #Construct data.
    Sig.probes <- read.csv(sprintf("%s/getMethdiff.%s.probes.significant.csv",
                                   dir.out,diff.dir),
                           stringsAsFactors=FALSE)[,1]
    ## Get nearby genes-----------------------
    nearGenes.file <- args[names(args) %in% "nearGenes"]
    if(length(nearGenes.file)==0){
      nearGenes.file <- sprintf("%s/%s.probes_nearGenes.rda",dir.out,diff.dir)
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
    params <- args[names(args) %in% c("percentage","permu.size","Pe","diffExp")]
    SigPair <- do.call(get.pair,c(list(mee=mee,probes=Sig.probes,
                                       nearGenes=nearGenes.file,
                                       permu.dir=permu.dir,
                                       dir.out=dir.out,
                                       cores=cores,label=diff.dir),
                                  params))
    
    ## promoter methylation correlation.
    # get promoter 
    promoter.probe <- suppressWarnings(get.feature.probe(promoter=TRUE, 
                                                         TSS.range=list(upstream=100,
                                                                        downstream=700)))
    mee <- fetch.mee(meth = meth.file, 
                     exp=exp.file, 
                     probeInfo=promoter.probe, 
                     geneInfo = geneAnnot,
                     TCGA=TRUE, 
                     genes=unique(SigPair$GeneID))
    params <- args[names(args) %in% "percentage"]
    Promoter.meth <- do.call(promoterMeth, c(list(mee=mee, sig.pvalue=0.01, save=FALSE),
                                             params))
    add <- SigPair[match(SigPair$GeneID, Promoter.meth$GeneID),"Raw.p"]
    SigPair <- cbind(SigPair, GSbPM.pvalue=add)
    write.csv(SigPair, file=sprintf("%s/getPair.%s.pairs.significant.csv",
                                    dir.out, diff.dir),
              row.names=FALSE)
    
    
    if(length(analysis) == 1){
      return(SigPair)
    } else {
      analysis <- analysis[!analysis %in% "pair"]
      rm(mee,SigPair)
      invisible(gc())
    }
  }
  
  # search enriched motif
  if("motif" %in% tolower(analysis)){
    message("###################\nMotif search\n###################\n\n")
    message(sprintf("Identify enriched motif for %smethylated probes",diff.dir))
    if(file.exists(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, diff.dir))){
      Sig.probes <- unique(read.csv(sprintf("%s/getPair.%s.pairs.significant.csv",
                                            dir.out, diff.dir), 
                                    stringsAsFactors=FALSE)$Probe)
    }else{
      stop(sprintf("%s/%s.pairs.significant.csv file doesn't exist",dir.out, diff.dir))
    }
    params <- args[names(args) %in% c("background.probes","lower.OR","min.incidence")]
    enriched.motif <- do.call(get.enriched.motif, c(list(probes=Sig.probes,
                                                         dir.out=dir.out,
                                                         label=diff.dir),
                                                    params))
    
    if(length(analysis) == 1){
      return(enriched.motif)
    } else {
      analysis <- analysis[!analysis %in% "motif"]
      rm(enriched.motif)
      invisible(gc())
    }
  }
  
  #search responsible TFs
  if("TF.search" %in% tolower(analysis)){
    message("###################\nSearch responsible TFs\n###################\n\n")
    ## make mee 
    #construct RNA seq data
    meth.file <- sprintf("%s/%s_meth_refined.rda",dir.out,disease)
    exp.file <- sprintf("%s/%s_RNA_refined.rda",dir.out,disease)
    probeInfo <- sprintf("%s/probeInfo_feature_distal.rda",dir.out)
    geneAnnot <- sprintf("%s/geneInfo.rda",dir.out)
    mee <- fetch.mee(meth=meth.file,exp=exp.file,TCGA=TRUE,
                     probeInfo=probeInfo,geneInfo=geneAnnot)
    message(sprintf("Identify regulatory TF for enriched motif in %smethylated probes",
                    diff.dir))
    enriched.motif <- args[names(args) %in% "enriched.motif"]
    if(length(enriched.motif)==0){
      enriched.motif <- sprintf("%s/getMotif.%s.enriched.motifs.rda",dir.out, diff.dir)
    }
    params <- args[names(args) %in% c("TFs", "motif.relavent.TFs","percentage")]
    TFs <- do.call(get.TFs, c(list(mee=mee, enriched.motif=enriched.motif,
                                   dir.out=dir.out, cores=cores, 
                                   label=diff.dir),
                              params))
  }
}

