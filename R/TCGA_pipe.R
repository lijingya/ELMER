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
#' @param group.col A column defining the groups of the sample. You can view the 
#' available columns using: colnames(MultiAssayExperiment::colData(data)).
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param wd A path shows working dirctory. Default is "./"
#' @param genes List of genes for which mutations will be verified. 
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
#' @importFrom SummarizedExperiment colData<-
#' @examples
#' \dontrun{
#'   distal.probe <- TCGA.pipe(disease = "LUSC", analysis="distal.enhancer", wd="~/")
#'   TCGA.pipe(disease = "LUSC",analysis = "all", genome = "hg19", cores = 1, permu.size=300, Pe=0.01)
#'   projects <- TCGAbiolinks:::getGDCprojects()$project_id
#'   projects <- gsub("TCGA-","",projects[grepl('^TCGA',projects,perl=TRUE)])
#'   for(proj in projects) TCGA.pipe(disease = proj,analysis = "download")
#'   plyr::alply(sort(projects),1,function(proj) {
#'        tryCatch({
#'          print(proj);
#'          TCGA.pipe(disease = proj,analysis = c("createMAE"))})
#'        }, .progress = "text")
#'   plyr::alply(sort(projects),1,function(proj) {
#'     tryCatch({
#'       print(proj);
#'       TCGA.pipe(disease = proj,
#'                  analysis = c("diffMeth","pair", "motif","TF.search"))})
#'   }, .progress = "text")
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
                                                        "Missense_Mutation",
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
  if (missing(disease)) 
    stop("Disease should be specified.\nDisease short name (such as LAML) 
         please check https://gdc-portal.nci.nih.gov")
  
  available.analysis <- c("download","distal.probes",
                          "createMAE","diffMeth","pair",
                          "motif","TF.search","all")
  # Replace all by all the other values
  if("all" %in% analysis[1] ) analysis <- grep("all", available.analysis, value = TRUE, invert = TRUE)
  
  if(any(!tolower(analysis) %in% tolower(analysis))) 
    stop(paste0("Availbale options for analysis argument are: ",
                paste(c("",available.analysis), collapse = "\n=> ")))
  
  disease <- toupper(disease)
  
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
  
  #---------------------------------------------------------------
  dir.out.root <- sprintf("%s/Result/%s",wd,disease)
  if(!file.exists(dir.out.root)) dir.create(dir.out.root, recursive = TRUE, showWarnings = FALSE)
  args <- list(...)
  
  if(any(sapply(analysis, function(x) tolower(x) %in% tolower(c("diffMeth","pair","motif","TF.search"))))) {
    if(!mode  %in% c("supervised","unsupervised")){
      stop("Set mode arugment to supervised or unsupervised")
    }
    if(mode %in% c("supervised")) {
      minSubgroupFrac <- 1
      message("=> ", mode, " was selected: using all samples")
    } else {
      minSubgroupFrac <- 0.2
      message("=> ", mode, " was selected: using ", minSubgroupFrac, " samples")
    }
    dir.out <- sprintf("%s/Result/%s/%s_%s_vs_%s/%s",wd,disease,group.col,group1,group2,diff.dir)
    message("=> Analysis results wil be save in:  ", dir.out)
    if(!file.exists(dir.out)) dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)
  }
  #-----------------------------------------------------
  ## select distal enhancer probes
  if(tolower("distal.probes") %in% tolower(analysis)){
    print.header("Select distal probes")
    params <- args[names(args) %in% c("TSS", "TSS.range","rm.chr")]
    params <- c(params,list("genome" = genome, "feature"= NULL))
    probeInfo <- do.call(get.feature.probe,params)
    save(probeInfo,file = sprintf("%s/probeInfo_distal_%s.rda",dir.out.root,genome))
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
      distal.probe <- sprintf("%s/probeInfo_distal_%s.rda",dir.out.root,genome)
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
      
      # if user set genes argument label Mutant WT will be added to mae
      if(!is.null(genes)) mae <- addMutCol(mae, disease, genes, mutant_variant_classification)
    } else {
      message("File already exists: ", file)
      mae <- get(load(file))
      if(!is.null(genes))  mae <- addMutCol(mae, disease, genes, mutant_variant_classification)
    }
    save(mae,file = file)
    message("File saved as: ", file)
    readr::write_tsv(as.data.frame(colData(mae)), path = sprintf("%s/%s_samples_info_%s.tsv",dir.out.root,disease,genome))
  }
  
  # Creates a record of the analysis and arguments called
  if(any(tolower(c("diffMeth","pair",
               "motif","TF.search")) %in% tolower(analysis))){
    createSummaryDocument(analysis = analysis, 
                          argument.values = args,
                          mae.path = sprintf("%s/%s_mae_%s.rda",dir.out.root,disease,genome),
                          genome = genome,
                          direction = diff.dir,
                          group.col = group.col,
                          group1 = group1,
                          group2 = group2,
                          results.path = dir.out)
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
    params <- c(params,list(diff.dir = diff.dir, 
                            dir.out = dir.out, 
                            cores = cores, 
                            minSubgroupFrac = minSubgroupFrac))
    diff.meth <- tryCatch({
      diff.meth <- do.call(get.diff.meth,c(params,list(data = mae,
                                                       group.col = group.col, 
                                                       group1 = group1,
                                                       group2 = group2)))
      diff.meth
    }, error = function(e) {
      message(e)
      return(NULL)
    })
    if(is.null(diff.meth)) return(NULL)
    
    if(length(analysis)==1) return(diff.meth)
  }
  
  # predict pair
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
    if(length(Sig.probes) == 0) {
      message("No significant probes were found")
      return(NULL)
    }
    
    message("Get nearby genes")
    file <- sprintf("%s/getPair.%s.pairs.significant.csv", dir.out, diff.dir)
    
    nearGenes.file <- args[names(args) %in% "nearGenes"]
    if(length(nearGenes.file)==0){
      nearGenes.file <- sprintf("%s/%s.probes_nearGenes.rda",dir.out,diff.dir)
      params <- args[names(args) %in% c("geneNum")]
      nearGenes <- do.call(GetNearGenes,
                           c(list(data = mae, 
                                  probes = Sig.probes,
                                  cores = cores),
                             params))
      save(nearGenes,file=nearGenes.file)
      message("File saved: ", nearGenes.file)
    } else {
      nearGenes.file <- nearGenes.file[["nearGenes"]]
    }
    
    # calculation
    message(sprintf("Identify putative probe-gene pair for %smethylated probes",diff.dir))
    
    # get pair
    permu.dir <- paste0(dir.out,"/permu")
    params <- args[names(args) %in% c("percentage","permu.size","Pe","raw.pvalue","diffExp","group.col")]
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
    
    # message("==== Promoter analysis ====")
    # message("calculate associations of gene expression with DNA methylation at promoter regions")
    # message("Fetching promoter regions")
    # file <- sprintf("%s/%s_mae_promoter_%s.rda",dir.out.root,disease, genome)
    # 
    # if(!file.exists(file)) {    
    #   ## promoter methylation correlation.
    #   # get promoter 
    #   suppressWarnings({
    #     promoter.probe <- get.feature.probe(promoter=TRUE, genome = genome,
    #                                         TSS.range=list(upstream = 200, downstream = 2000))
    #   })
    #   
    #   if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    #   meth.file <- sprintf("%s/%s_meth_%s.rda",Data,disease, genome)
    #   if(is.null(Data)) Data <- sprintf("%s/Data/%s",wd,disease)
    #   exp.file <- sprintf("%s/%s_RNA_%s.rda",Data,disease, genome)
    #   
    #   mae.promoter <- createMAE(met           = meth.file, 
    #                             exp           = exp.file, 
    #                             filter.probes = promoter.probe,
    #                             genome        = genome,
    #                             met.platform  = "450K",
    #                             linearize.exp = TRUE,
    #                             save          = FALSE,
    #                             TCGA          = TRUE)
    #   if(!all(sample.type %in% colData(mae)[,group.col])){
    #     message("There are no samples for both groups")
    #     return(NULL)
    #   }
    #   mae.promoter <- mae.promoter[,colData(mae.promoter)[,group.col] %in% sample.type]
    #   save(mae.promoter,file = file)
    # } else {
    #   mae.promoter <- get(load(file))
    # }
    # params <- args[names(args) %in% "percentage"]
    # Promoter.meth <- do.call(promoterMeth, c(list(data=mae.promoter, sig.pvalue=0.01, save=FALSE),
    #                                          params))
    # write.csv(Promoter.meth, 
    #           file = sprintf("%s/promoter.%s.analysis.csv", dir.out, diff.dir),
    #           row.names=FALSE)
    # add <- SigPair[match(SigPair$GeneID, Promoter.meth$GeneID),"Raw.p"]
    # SigPair <- cbind(SigPair, GSbPM.pvalue = add)
    if(is.null(SigPair)) {
      message("No significant pair probe genes found")
      return(NULL)
    }
    write.csv(SigPair, 
              file = sprintf("%s/getPair.%s.pairs.significant.csv", dir.out, diff.dir),
              row.names=FALSE)
    
    if(length(analysis) == 1) return(SigPair)
  }
  
  # search enriched motif
  if("motif" %in% tolower(analysis)){
    print.header("Motif search")
    mae.file <- sprintf("%s/%s_mae_%s.rda",dir.out.root,disease,genome)
    if(!file.exists(mae.file)){
      message("MAE not found, please run pipe with createMAE or all options")
      return(NULL)
    }
    load(mae.file)
    message(sprintf("Identify enriched motif for %smethylated probes",diff.dir))
    if(file.exists(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, diff.dir))){
      Sig.probes <- readr::read_csv(sprintf("%s/getPair.%s.pairs.significant.csv", dir.out, diff.dir), 
                                    col_names = TRUE,
                                    col_types = c("cccicdd"))
      if(length(unique(Sig.probes)) < 10) {
        message ("No significants pairs were found in the previous step") 
        return(NULL)
      }
      Sig.probes <- unique(Sig.probes$Probe)
      
    } else {
      message(sprintf("%s/%s.pairs.significant.csv file doesn't exist",dir.out, diff.dir))
      return(NULL)
    }
    params <- args[names(args) %in% c("background.probes","lower.OR","min.incidence")]
    
    newenv <- new.env()
    if(genome == "hg19") data("Probes.motif.hg19.450K", package = "ELMER.data", envir = newenv)
    if(genome == "hg38") data("Probes.motif.hg38.450K", package = "ELMER.data", envir = newenv)
    probes.motif <- get(ls(newenv)[1],envir=newenv)   
    
    enriched.motif <- do.call(get.enriched.motif, 
                              c(list(data         = mae,
                                     probes.motif = probes.motif,
                                     probes       = Sig.probes,
                                     dir.out      = dir.out,
                                     label        = diff.dir,
                                     plot.title   = paste0("OR for paired probes ",
                                                           diff.dir, " methylated in ",
                                                           group1, " vs ",group2, "(group: ",group.col,")")),
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

#' Adds mutation information to MAE
#' @param data MAE object
#' @param disease TCGA disease (LUSC, GBM, etc)
#' @param genes list of genes to add information
#' @param mutant_variant_classification List of mutant_variant_classification that will be 
#' consider a sample mutant or not.
#' @examples
#' \dontrun{
#'  data <- ELMER:::getdata("elmer.data.example") # Get data from ELMER.data
#'  data <- addMutCol(data, "LUSC","TP53")
#' }
addMutCol <- function(data, 
                      disease, 
                      genes, 
                      mutant_variant_classification = c("Frame_Shift_Del",
                                                        "Frame_Shift_Ins",
                                                        "Missense_Mutation",
                                                        "Nonsense_Mutation",
                                                        "Splice_Site",
                                                        "In_Frame_Del",
                                                        "In_Frame_Ins",
                                                        "Translation_Start_Site",
                                                        "Nonstop_Mutation")){
  maf <- TCGAbiolinks::GDCquery_Maf(disease , pipeline = "mutect2")
  for(gene in genes) {
    if(gene %in% maf$Hugo_Symbol) {
      message("Adding information for gene: ", gene)
      aux <- maf %>% filter(Hugo_Symbol == gene) # Select only mutation on that gene
      idx <- unique(unlist(sapply(mutant_variant_classification,function(x) grep(x,aux$Variant_Classification, ignore.case = TRUE))))
      aux <- aux[idx,]
      mutant.samples <- substr(aux$Tumor_Sample_Barcode,1,16)
      colData(data)[,gene] <- "Normal"
      colData(data)[colData(data)$TN == "Tumor", gene] <- "WT"
      colData(data)[colData(data)$TN == "Tumor" & 
                      colData(data)$sample %in% mutant.samples,gene] <- "Mutant"
      message("The column ", gene, " was create in the MAE object")
      print(plyr::count(colData(data)[,gene]))
    } else {
      message("No mutation found for: ", gene)
    }
  }
  return(data)
}

#' @title Create summary document for TCGA.pipe function
#' @description This function will create a text file with the 
#' date of the last run, which aanalysis were performed, the values of
#' the arguments so the user can keep track 
createSummaryDocument <- function(analysis = "all", 
                                  argument.values = "defaults",
                                  genome = NULL,
                                  mae.path = NULL,
                                  direction = NULL,
                                  group.col = NULL,
                                  group1 = NULL,
                                  group2 = NULL,
                                  results.path = NULL){
  print("Recording analysis information into TCGA.pipe_records.txt")
  df <- paste0("oooooooooooooooooooooooooooooooooooo\n",
               "o date: ",Sys.time(),"\n",
               "o analysis: ", paste(analysis, collapse = ","), "\n",
               "o genome: ",  ifelse(is.null(genome),"",genome), "\n",
               "o mae.path: ",  ifelse(is.null(mae.path),"",mae.path), "\n",
               "o direction: ", ifelse(is.null(direction),"",direction), "\n",
               "o group.col: ",  ifelse(is.null(group.col),"",group.col), "\n",
               "o group1: ",   ifelse(is.null(group1),"",group1), "\n",
               "o group2: ",   ifelse(is.null(group2),"",group2), "\n",
               "o results.path: ",  ifelse(is.null(results.path),"",results.path), "\n",
               "o argument.values: ",paste(paste0(names(argument.values),"=",as.character(argument.values)), collapse = ",")
  )
  fileConn <- file(file.path(results.path,"TCGA.pipe_records.txt"),open = "a+")
  write(df, fileConn, append=TRUE)
  close(fileConn)
}