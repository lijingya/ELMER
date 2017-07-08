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
#' @param gene List of genes for which mutations will be verified. 
#' A column in the MAE with the name of the gene
#' will be created with two groups WT (tumor samples without mutation), MUT (tumor samples w/ mutation), 
#' NA (not tumor samples)
#' @param mode This option will automatically set the percentage of samples to be used in the analysis.
#' Options: "supervised" (use 100\% of samples) or "unsupervised" (use 20\% of samples).
#' @param cores A interger which defines number of core to be used in parallel process. 
#' Default is 1: don't use parallel process.
#' @param diff.dir A character can be "hypo" or "hyper", showing dirction DNA methylation changes.  
#' If it is "hypo", get.diff.meth function will identify all significantly hypomethylated
#' CpG sites; If "hyper", get.diff.meth function will identify all significantly hypermethylated
#' CpG sites
#' @param mutant_variant_classification List of TCGA variant classification from MAF files to consider a samples
#' mutant. Only used when argument gene is set.
#' @param Data A path shows the folder containing DNA methylation, expression and clinic data
#' @param ... A list of parameters for functions: GetNearGenes, get.feature.probe, 
#' get.diff.meth, get.pair
#' @return Different analysis results.
#' @export 
#' @examples
#' \dontrun{
#'   distal.probe <- TCGA.pipe(disease = "LUSC", analysis="distal.enhancer", wd="~/")
#'   TCGA.pipe(disease = "LUSC",analysis = "all", genome = "hg19", cores = 1, permu.size=300, Pe=0.01)
#'   projects <- TCGAbiolinks:::getGDCprojects()$project_id
#'   projects <- gsub("TCGA-","",projects[grepl('^TCGA',projects,perl=TRUE)])
#'   for(proj in projects) TCGA.pipe(disease = proj,analysis = "download")
#'   plyr::alply(sort(projects),1,function(proj) {tryCatch({print(proj);TCGA.pipe(disease = proj,analysis = c("createMAE"))})}, .progress = "text")
#'   plyr::alply(sort(projects),1,function(proj) {tryCatch({print(proj);TCGA.pipe(disease = proj,analysis = c("diffMeth","pair", "motif","TF.search"))})}, .progress = "text")
#'
#'   # Evaluation mutation
#'   TCGA.pipe(disease = "LUSC",analysis = "createMAE",gene = "NFE2L2")
#'   TCGA.pipe(disease = "LUSC",analysis = c("diffMeth","pair", "motif","TF.search"), 
#'             mode = "supervised",
#'             group.col = "NFE2L2", group1 = "Mutant", group2 = "WT",
#'             diff.dir = c("hypo"),
#'             dir.out = "LUSC_NFE2L2_MutvsWT")
#' }
TCGA.pipe <- function(disease,
                      genome = "hg38",
                      analysis = "all",
                      wd = "./",
                      cores = 1,
                      mode = "unsupervised",
                      Data = NULL, 
                      diff.dir = "hypo",
                      genes = NULL,
                      mutant_variant_classification = c("Frame_Shift_Del",
                                                 "Frame_Shift_Ins",
                                                 " Missense_Mutation",
                                                 "Nonsense_Mutation",
                                                 "Splice_Site",
                                                 "In_Frame_Del",
                                                 "In_Frame_Ins",
                                                 "Translation_Start_Site",
                                                 "Nonstop_Mutation"),
                      group.col = "TN", 
                      group1 = "Tumor",
                      group2 = "Normal",
                      ...){
  if(!mode  %in% c("supervised","unsupervised")){
    stop("Set mode arugment to supervised or unsupervised")
  }
  if(mode %in% c("supervised")) {
    minSubgroupFrac <- 1
    message("=> ",mode, " was selected: using all samples")
  } else {
    minSubgroupFrac <- 0.2
    message("=> ",mode, " was selected: using ", minSubgroupFrac, " samples")
  }
  if (missing(disease)) 
    stop("Disease should be specified.\nDisease short name (such as LAML) 
         please check https://gdc-portal.nci.nih.gov")
  
  available.analysis <- c("download","distal.probes",
                          "createMAE","diffMeth","pair",
                          "motif","TF.search","all")
  if (analysis[1] == "all") analysis <- grep("all", available.analysis, value = TRUE, invert = TRUE)
  
  if(any(!tolower(analysis) %in% tolower(analysis))) 
    stop(paste0("Availbale options for analysis argument are: ",
                paste(c("",available.analysis), collapse = "\n=> ")))
  
  disease <- toupper(disease)
  dir.out.root <- sprintf("%s/Result/%s",wd,disease)
  dir.out <- sprintf("%s/Result/%s/%s_%s_vs_%s/%s",wd,disease,group.col,group1,group2,diff.dir)
  message("=> Saving results to ", dir.out)
  if(!file.exists(dir.out)) dir.create(dir.out,recursive = TRUE, showWarnings = FALSE)
  args <- list(...)
  
  # Download 
  if("download" %in% tolower(analysis)){
    print.header("Download data")
    if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    params <- c()
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
    file <- sprintf("%s/%s_mae_%s.rda",dir.out.root,disease,genome)
    if(!file.exists(file)) {
      sample.type <- c(group1,group2)
      
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
                       save          = FALSE,
                       linearize.exp = TRUE,
                       TCGA          = TRUE)
      
      # if user set genes argument label MUT WT will be added to mae
      if(!is.null(genes)) {
        maf <- TCGAbiolinks::GDCquery_Maf(disease , pipeline = "mutect2")
        for(g in genes) {
          if(g %in% maf$Hugo_Symbol) {
            message("Adding information for gene: ", g)
            aux <- filter(maf, Hugo_Symbol == g & Variant_Classification %in% mutant_variant_classification)
            mutant.samples <- substr(aux$Tumor_Sample_Barcode,1,15)
            colData(mae)[,g] <- NA
            colData(mae)[colData(mae)$TN == " Tumor",g] <- "WT"
            colData(mae)[colData(mae)$TN == " Tumor" & 
                           colData(mae)$samples %in% mutant.samples,g] <- paste0(g, " mutant")
            print(plyr::count(colData(mae)[,g]))
          } else {
            message("No mutation found for: ", g)
          }
        }
      }
      
      if(!all(sample.type %in% colData(mae)[,group.col])){
        print(table(colData(mae)[,group.col]))
        print(mae)
        message("There are no samples for both groups")
        return(NULL)
      }
      
      mae <- mae[,colData(mae)[,group.col] %in% sample.type]
      save(mae,file = file)
      message("File saved as: ", file)
      readr::write_tsv(as.data.frame(colData(mae)), path = sprintf("%s/%s_samples_info_%s.tsv",dir.out,disease,genome))
    } else {
      message("File already exists: ", file)
      if(!is.null(genes)) {
        mae <- get(load(file))
        maf <- TCGAbiolinks::GDCquery_Maf(disease , pipeline = "mutect2")
        for(g in genes) {
          if(g %in% maf$Hugo_Symbol) {
            message("Adding information for gene: ", g)
            aux <- filter(maf, Hugo_Symbol == g & !Variant_Classification %in% c("3'UTR","3'Flank","5'UTR","5'Flank","Silent","Intron"))
            mutant.samples <- substr(aux$Tumor_Sample_Barcode,1,16)
            colData(mae)[,g] <- NA
            colData(mae)[colData(mae)$TN == "Tumor",g] <- "WT"
            colData(mae)[colData(mae)$TN == "Tumor" & 
                           colData(mae)$sample %in% mutant.samples,g] <- paste0("Mutant")
            print(plyr::count(colData(mae)[,g]))
          } else {
            message("No mutation found for: ", g)
          }
        }
        save(mae,file = file)
      }
    }
  }
  
  # get differential DNA methylation
  if(tolower("diffMeth") %in% tolower(analysis)){
    print.header("Get differential DNA methylation loci")
    mae.file <- sprintf("%s/%s_mae_%s.rda",dir.out.root,disease,genome)
    if(!file.exists(mae.file)){
      message("MAE not found, please run pipe with createMAE or all options")
      return(NULL)
    }
    load(mae.file)
    params <- args[names(args) %in% c("pvalue","sig.dif")]
    params <- c(params,list(diff.dir=diff.dir, dir.out=dir.out, cores=cores, minSubgroupFrac = minSubgroupFrac))
    diff.meth <- do.call(get.diff.meth,c(params,list(data = mae,
                                                     group.col = group.col, 
                                                     group1 = group1,
                                                     group2 = group2)))
    if(length(analysis)==1) return(diff.meth)
  }
  
  #predict pair
  if("pair" %in% tolower(analysis)){
    print.header("Predict pairs")
    mae.file <- sprintf("%s/%s_mae_%s.rda",dir.out.root,disease,genome)
    if(!file.exists(mae.file)){
      message("MAE not found, please run pipe with createMAE or all options")
      return(NULL)
    }
    load(mae.file)
    
    Sig.probes <- read.csv(sprintf("%s/getMethdiff.%s.probes.significant.csv",
                                   dir.out,diff.dir),
                           stringsAsFactors=FALSE)[,1]
    if(length(Sig.probes) == 0) stop("No significant probes were found")
    ## Get nearby genes-----------------------
    message("Get nearby genes")
    file <- sprintf("%s/getPair.%s.pairs.significant.csv", dir.out, diff.dir)
    
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
      params <- args[names(args) %in% c("percentage","permu.size","Pe","pvalue","diffExp","group.col")]
      SigPair <- do.call(get.pair,
                         c(list(data      = mae,
                                nearGenes = nearGenes.file,
                                permu.dir = permu.dir,
                                group.col = group.col, 
                                group1    = group1,
                                minSubgroupFrac = min(1,minSubgroupFrac * 2),
                                group2    = group2,
                                dir.out   = dir.out,
                                cores     = cores,
                                label     = diff.dir),
                           params))
      
    message("==== Promoter analysis ====")
    message("calculate associations of gene expression with DNA methylation at promoter regions")
    message("Fetching promoter regions")
    file <- sprintf("%s/%s_mae_promoter_%s.rda",dir.out.root,disease, genome)
    
    if(!file.exists(file)){    
      ## promoter methylation correlation.
      # get promoter 
      suppressWarnings({
        promoter.probe <- get.feature.probe(promoter=TRUE, genome = genome,
                                            TSS.range=list(upstream = 200, downstream = 2000))
      })
     
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
                                save          = FALSE,
                                TCGA          = TRUE)
      if(!all(sample.type %in% colData(mae)[,group.col])){
        message("There are no samples for both groups")
        return(NULL)
      }
      mae.promoter <- mae.promoter[,colData(mae.promoter)[,group.col] %in% sample.type]
      save(mae.promoter,file = file)
    } else {
      mae.promoter <- get(load(file))
    }
    params <- args[names(args) %in% "percentage"]
    Promoter.meth <- do.call(promoterMeth, c(list(data=mae.promoter, sig.pvalue=0.01, save=FALSE),
                                             params))
    write.csv(Promoter.meth, 
              file = sprintf("%s/promoter.%s.analysis.csv", dir.out, diff.dir),
              row.names=FALSE)
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
    mae.file <- sprintf("%s/%s_mae_%s.rda",dir.out.root,disease,genome)
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
                          group.col      = group.col, 
                          group1         = group1,
                          group2         = group2,
                          minSubgroupFrac =  min(1,minSubgroupFrac * 2),
                          enriched.motif = enriched.motif,
                          dir.out        = dir.out, 
                          cores          = cores, 
                          label          = diff.dir),
                     params))
    if(length(analysis) == 1) return(TFs)
  }
  
}

