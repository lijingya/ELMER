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
#'   distal.probe <- TCGA.pipe(disease = "LUSC", analysis="distal.enhancer", wd="~/")
#'   TCGA.pipe(disease = "LUSC",analysis = "all", genome = "hg19", cores = 1, permu.size=300, Pe=0.01)
#' }
TCGA.pipe <- function(disease,
                      genome = "hg38",
                      analysis = "all",
                      wd = "./",
                      cores = 1,
                      Data = NULL, 
                      diff.dir = "hypo",
                      ...){
  if (missing(disease)) 
    stop("Disease should be specified.\nDisease short name (such as LAML) 
         please check https://gdc-portal.nci.nih.gov")
  
  available.analysis <- c("download","distal.probes",
                          "createMAE","diffMeth","pair",
                          "motif","TF.search","all")
  if (analysis[1] == "all") analysis <- grep("all", available.analysis, value = T, invert = T)
  
  if(any(!tolower(analysis) %in% tolower(analysis))) 
    stop(paste0("Availbale options for analysis argument are: ",
                paste(c("",available.analysis), collapse = "\n=> ")))
  
  disease <- toupper(disease)
  dir.out <- sprintf("%s/Result/%s",wd,disease)
  if(!file.exists(dir.out)) dir.create(dir.out,recursive = TRUE)
  args <- list(...)
  
  # Download 
  if("download" %in% tolower(analysis)){
    print.header("Download data")
    if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    params <- args[names(args) %in% c("RNAtype","Methfilter")]
    params$disease <- disease
    params$basedir <- sprintf("%s/Data",wd)
    params$genome <- genome
    do.call(getTCGA,params)
    analysis <- analysis[!analysis %in% "download"]
  }
  
  ## select distal enhancer probes
  if(tolower("distal.probes") %in% tolower(analysis)){
    print.header("Select distal probes")
    params <- args[names(args) %in% c("TSS", "TSS.range","rm.chr")]
    params <- c(params,list("genome" = genome, "feature"= NULL))
    probeInfo <- do.call(get.feature.probe,params)
    save(probeInfo,file = sprintf("%s/probeInfo_distal_%s.rda",dir.out,genome))
    if(length(analysis) == 1){
      return(probeInfo)
    } else {
      analysis <- analysis[!analysis %in% "distal.probes"]
      invisible(gc())
    }
  }
  
  if(tolower("createMAE") %in% tolower(analysis)){
    print.header("Creating Multi Assay Experiment")

    group.col <- "TN"
    sample.type <- c("Tumor","Normal")
    
    if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    meth.file <- sprintf("%s/%s_meth_%s.rda",Data,disease,genome)
    if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    exp.file <- sprintf("%s/%s_RNA_%s.rda",Data,disease,genome)

    ## get distal probe info
    distal.probe <- sprintf("%s/probeInfo_distal_%s.rda",dir.out,genome)
    if(!file.exists(distal.probe)){
      params <- args[names(args) %in% c("TSS","TSS.range","rm.chr")]
      params <- c(params,list("genome" = genome, "feature"= NULL))
      distal.probe <- suppressWarnings(do.call(get.feature.probe,params))
    }
    
    mae <- createMAE(met           = meth.file, 
                     exp           = exp.file, 
                     filter.probes = distal.probe,
                     genome        = genome,
                     met.platform  = "450K",
                     save = FALSE,
                     linearize.exp = TRUE,
                     TCGA          = TRUE)
    if(!all(sample.type %in% colData(mae)[,group.col])){
      message("There are no samples for both groups")
      return(NULL)
    }
    mae <- mae[,colData(mae)[,group.col] %in% sample.type]
    save(mae,file = sprintf("%s/%s_mae_%s.rda",dir.out,disease,genome))
    readr::write_tsv(as.data.frame(colData(mae)), path = sprintf("%s/%s_samples_info_%s.tsv",dir.out,disease,genome))
  }
  
  # get differential DNA methylation
  if(tolower("diffMeth") %in% tolower(analysis)){
    print.header("Get differential DNA methylation loci")
    mae.file <- sprintf("%s/%s_mae_%s.rda",dir.out,disease,genome)
    if(!file.exists(mae.file)){
      message("MAE not found, please run pipe with createMAE or all options")
      return(NULL)
    }
    load(mae.file)
    params <- args[names(args) %in% c("percentage","pvalue","sig.dif")]
    params <- c(params,list(diff.dir=diff.dir, dir.out=dir.out, cores=cores))
    diff.meth <- do.call(get.diff.meth,c(params,list(data=mae,group.col = "TN")))
    if(length(analysis)==1) return(diff.meth)
  }
  
  #predict pair
  if("pair" %in% tolower(analysis)){
    print.header("Predict pairs")
    mae.file <- sprintf("%s/%s_mae_%s.rda",dir.out,disease, genome)
    if(!file.exists(mae.file)){
      message("MAE not found, please run pipe with createMAE or all options")
      return(NULL)
    }
    load(mae.file)
    
    Sig.probes <- read.csv(sprintf("%s/getMethdiff.%s.probes.significant.csv",
                                   dir.out,diff.dir),
                           stringsAsFactors=FALSE)[,1]
    ## Get nearby genes-----------------------
    message("Get nearby genes")
    
    nearGenes.file <- args[names(args) %in% "nearGenes"]
    if(length(nearGenes.file)==0){
      nearGenes.file <- sprintf("%s/%s.probes_nearGenes.rda",dir.out,diff.dir)
      if(!file.exists(nearGenes.file)){
        params <- args[names(args) %in% c("geneNum")]
        nearGenes <- do.call(GetNearGenes,
                             c(list(data = mae, 
                                    probes = Sig.probes,
                                    cores = cores),
                               params))
        save(nearGenes,file=nearGenes.file)
      }
    } else {
      nearGenes.file <- nearGenes.file[["nearGenes"]]
    }
    ## calculation
    message(sprintf("Identify putative probe-gene pair for %smethylated probes",diff.dir))
    
    ## get pair
    permu.dir <- paste0(dir.out,"/permu")
    params <- args[names(args) %in% c("percentage","permu.size","Pe","diffExp","calculate.Pe","group.col")]
    SigPair <- do.call(get.pair,
                       c(list(data      = mae,
                              nearGenes = nearGenes.file,
                              permu.dir = permu.dir,
                              group.col = "TN",
                              dir.out   = dir.out,
                              cores     = cores,
                              label     = diff.dir),
                         params))
    
    message("calculate associations of gene expression with DNA methylation at promoter regions")
    message("Fetching promoter regions")
    
    ## promoter methylation correlation.
    # get promoter 
    suppressWarnings({
      promoter.probe <- get.feature.probe(promoter=TRUE, 
                                          TSS.range=list(upstream=100, downstream=700))
    })
    group.col <- "TN"
    sample.type <- c("Tumor","Normal")
    
    if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    meth.file <- sprintf("%s/%s_meth_%s.rda",Data,disease, genome)
    if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    exp.file <- sprintf("%s/%s_RNA_%s.rda",Data,disease, genome)
    
    mae.promoter <- createMAE(met           = meth.file, 
                              exp           = exp.file, 
                              filter.probes = promoter.probe,
                              genome        = genome,
                              met.platform  = "450K",
                              linearize.exp = TRUE,
                              TCGA          = TRUE)
    if(!all(sample.type %in% colData(mae)[,group.col])){
      message("There are no samples for both groups")
      return(NULL)
    }
    mae.promoter <- mae.promoter[,colData(mae.promoter)[,group.col] %in% sample.type]
    save(mae.promoter,file = sprintf("%s/%s_mae_promoter_%s.rda",dir.out,disease, genome))
    
    params <- args[names(args) %in% "percentage"]
    Promoter.meth <- do.call(promoterMeth, c(list(data=mae.promoter, sig.pvalue=0.01, save=FALSE),
                                             params))
    add <- SigPair[match(SigPair$GeneID, Promoter.meth$GeneID),"Raw.p"]
    SigPair <- cbind(SigPair, GSbPM.pvalue = add)
    write.csv(SigPair, 
              file = sprintf("%s/getPair.%s.pairs.significant.csv", dir.out, diff.dir),
              row.names=FALSE)
    
    if(length(analysis) == 1) return(SigPair)
  }
  
  # search enriched motif
  if("motif" %in% tolower(analysis)){
    print.header("Motif search")

    message(sprintf("Identify enriched motif for %smethylated probes",diff.dir))
    if(file.exists(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, diff.dir))){
      Sig.probes <- unique(read.csv(sprintf("%s/getPair.%s.pairs.significant.csv",
                                            dir.out, diff.dir), 
                                    stringsAsFactors=FALSE)$Probe)
    } else {
      stop(sprintf("%s/%s.pairs.significant.csv file doesn't exist",dir.out, diff.dir))
    }
    params <- args[names(args) %in% c("background.probes","lower.OR","min.incidence")]
    
    newenv <- new.env()
    if(genome == "hg19") data("Probes.motif.hg19.450K", package = "ELMER.data", envir = newenv)
    if(genome == "hg38") data("Probes.motif.hg38.450K", package = "ELMER.data", envir = newenv)
    probes.motif <- get(ls(newenv)[1],envir=newenv)   
    
    enriched.motif <- do.call(get.enriched.motif, 
                              c(list(probes.motif = probes.motif,
                                     probes       = Sig.probes,
                                     dir.out      = dir.out,
                                     label        = diff.dir),
                                params))
    
    if(length(analysis) == 1) return(enriched.motif)
  }
  
  #search responsible TFs
  if(tolower("TF.search") %in% tolower(analysis)){
    print.header("Search responsible TFs")
    ## load mae
    mae.file <- sprintf("%s/%s_mae.rda",dir.out,disease)
    if(!file.exists(mae.file)){
      message("MAE not found, please run pipe with createMAE or all options")
      return(NULL)
    }
    load(mae.file)
    #construct RNA seq data
    print.header(sprintf("Identify regulatory TF for enriched motif in %smethylated probes",
                    diff.dir), "subsection")
    enriched.motif <- args[names(args) %in% "enriched.motif"]
    if(length(enriched.motif) == 0){
      enriched.motif <- sprintf("%s/getMotif.%s.enriched.motifs.rda", dir.out, diff.dir)
    }
    params <- args[names(args) %in% c("TFs", "motif.relavent.TFs","percentage")]
    TFs <- do.call(get.TFs, 
                   c(list(data           = mae, 
                          enriched.motif = enriched.motif,
                          dir.out        = dir.out, 
                          cores          = cores, 
                          label          = diff.dir),
                     params))
    if(length(analysis) == 1) return(TFs)
  }
           
}

