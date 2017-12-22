#' @title get.feature.probe to select probes within promoter regions or distal regions.
#' @description 
#' get.feature.probe is a function to select the probes falling into 
#' distal feature regions or promoter regions.
#' @importFrom GenomicRanges promoters 
#' @description This function selects the probes on HM450K that either overlap 
#' distal biofeatures or TSS promoter. 
#' @param promoter A logical.If TRUE, function will ouput the promoter probes.
#' If FALSE, function will ouput the distal probes overlaping with features. The 
#' default is FALSE.
#' @param met.platform DNA methyaltion platform to retrieve data from: EPIC or 450K (default)
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @param feature A GRange object containing biofeature coordinate such as 
#' enhancer coordinates. 
#' If NULL only distal probes (2Kbp away from TSS will be selected)
#' feature option is only usable when promoter option is FALSE.
#' @param TSS A GRange object contains the transcription start sites. When promoter is FALSE, Union.TSS
#' in \pkg{ELMER.data} will be used for default. When promoter is TRUE, UCSC gene TSS will
#' be used as default (see detail). User can specify their own preference TSS annotation. 
#' @param TSS.range A list specify how to define promoter regions. 
#' Default is upstream =2000bp and downstream=2000bp.
#' @param rm.chr A vector of chromosome need to be remove from probes such as chrX chrY or chrM
#' @return A GRange object containing probes that satisfy selecting critiria.
#' @export 
#' @importFrom S4Vectors queryHits subjectHits
#' @details 
#'  In order to get real distal probes, we use more comprehensive annotated TSS by both 
#'  GENCODE and UCSC. However, to get probes within promoter regions need more
#'  accurate annotated TSS such as UCSC. Therefore, there are different settings for
#'  promoter and distal probe selection. But user can specify their own favorable
#'  TSS annotation. Then there won't be any difference between promoter and distal
#'  probe selection.
#'  @return A GRanges object contains the coordinate of probes which locate 
#'  within promoter regions or distal feature regions such as union enhancer from REMC and FANTOM5.
#'  @usage get.feature.probe(feature, 
#'                           TSS, 
#'                           TSS.range = list(upstream = 2000, downstream = 2000), 
#'                           promoter = FALSE, rm.chr = NULL)
#' @examples 
#' # get distal enhancer probe
#' \dontrun{
#' Probe <- get.feature.probe()
#' }
#' # get promoter probes
#' \dontrun{
#' Probe <- get.feature.probe(promoter=FALSE)
#' }
#' # get distal enhancer probe remove chrX chrY
#' Probe2 <- get.feature.probe(rm.chr=c("chrX", "chrY"))
get.feature.probe <- function(feature = NULL,
                              TSS,
                              genome = "hg38",
                              met.platform = "450K",
                              TSS.range = list(upstream = 2000, downstream = 2000),
                              promoter = FALSE,
                              rm.chr = NULL){
  probe <- getInfiniumAnnotation(toupper(met.platform),genome)
  # We will rmeove the rs probes, as they should not be used in the analysis
  probe <- probe[!grepl("rs",names(probe)),]
  probe <- probe[!probe$MASK.general,] # remove masked probes
  if(!is.null(rm.chr)) probe <- probe[!as.character(seqnames(probe)) %in% rm.chr]
  
  if(missing(TSS)){
    # The function getTSS gets the transcription coordinantes from Ensemble (GENCODE)
    TSS <- getTSS(genome = genome)
  }
  suppressWarnings({
    promoters <- promoters(TSS,
                           upstream = TSS.range[["upstream"]], 
                           downstream = TSS.range[["downstream"]])
  })
  
  if(!promoter){
    probe <- probe[setdiff(1:length(probe),unique(queryHits(findOverlaps(probe,promoters))))]
    
    
    if(is.null(feature)) {
      message("Returning distal probes: ", length(probe))
      return(probe)
    }
    if(is(feature,"GRanges")) {             
      probe <- probe[unique(queryHits(findOverlaps(probe,feature)))]
      message("Returning distal probes overlapping with features: ", length(probe))
      
    } else {
      stop("feature is not GRanges object.")
    }
  } else {
    probe <- probe[unique(queryHits(findOverlaps(probe,promoters)))]
  }
  return(probe)
}


## get differential methylated probes-------------------------
## TCGA pipe don't specify dir.out
#' get.diff.meth to identify hypo/hyper-methylated CpG sites on HM450K between control and experimental 
#' groups such as normal verus tumor samples.
#' @description 
#' get.diff.meth applys one-way t-test to identify the CpG sites that are significantly 
#' hypo/hyper-methyalated using proportional samples (defined by minSubgroupFrac option) from group 1 
#' and group 2. The P values will be adjusted by Benjamini-Hochberg method. 
#' Option pvalue and sig.dif will be the criteria (cutoff) for selecting significant 
#' differentially methylated CpG sites.
#' If save is TURE, two getMethdiff.XX.csv files will be generated (see detail).
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. 
#' See \code{\link{createMAE}} function.
#' @param group.col A column defining the groups of the sample. You can view the 
#' available columns using: colnames(MultiAssayExperiment::colData(data)).
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param diff.dir A character can be "hypo" or "hyper", showing differential 
#' methylation dirction.  It can be "hypo" which is only selecting hypomethylated probes; 
#' "hyper" which is only selecting hypermethylated probes; 
#' @param cores A interger which defines the number of cores to be used in parallel 
#' process. Default is 1: no parallel process.
#' @param minSubgroupFrac A number ranging from 0 to 1, 
#' specifying the fraction of extreme samples from group 1 and group 2 
#' that are used to identify the differential DNA methylation. 
#' The default is 0.2 because we typically want to be able to detect a specific 
#' (possibly unknown) molecular subtype among tumor; these subtypes often make up only 
#' a minority of samples, and 20\% was chosen as a lower bound for the purposes of statistical power. 
#' If you are using pre-defined group labels, such as treated replicates vs. untreated replicated, 
#' use a value of 1.0 (Supervised mode)
#' @param pvalue A number specifies the significant P value (adjusted P value by BH) 
#' threshold Limit for selecting significant hypo/hyper-methylated probes. Default is 0.01
#' If pvalue is smaller than pvalue than it is considered significant.
#' @param sig.dif A number specifies the smallest DNA methylation difference as a cutoff for 
#' selecting significant hypo/hyper-methylated probes. Default is 0.3.
#' @param dir.out A path specify the directory for outputs. Default is is current directory.
#' @param test Statistical test to be used. Options: t.test (DEFAULT), wilcox.test
#' @param save A logic. When TRUE, two getMethdiff.XX.csv files will be generated (see detail)
#' @param min.samples Minimun number of samples to use in the analysis. Default 5.
#' If you have 10 samples in one group, minSubgroupFrac is 0.2 this will give 2 samples 
#' in the lower quintile, but then 5 will be used.
#' @details 
#'  save: 
#'  When save is TRUE, function will generate two XX.csv files.The first one is named 
#'  getMethdiff.hypo.probes.csv (or getMethdiff.hyper.probes.csv depends on diff.dir). 
#'  The first file contains all statistic results for each probe. Based on this
#'  file, user can change different P value or sig.dir cutoff to select the significant results
#'  without redo the analysis. The second file is named getMethdiff.hypo.probes.significant.csv
#'  (or getMethdiff.hyper.probes.significant.csv depends on diff.dir). This file contains
#'  statistic results for the probes that pass the significant criteria (P value and sig.dir).
#'  When save is FALSE, a data frame R object will be generate which contains the same
#'  information with the second file.
#' @return Statistics for all probes and significant hypo or hyper-methylated probes.
#' @export 
#' @importFrom readr write_csv
#' @importFrom plyr adply
#' @importFrom stats p.adjust
#' @importFrom MultiAssayExperiment colData
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' data <- ELMER:::getdata("elmer.data.example")
#' Hypo.probe <- get.diff.meth(data, 
#'                             diff.dir="hypo",
#'                             group.col = "definition", 
#'                             group1 = "Primary solid Tumor", 
#'                             group2 = "Solid Tissue Normal",
#'                             sig.dif = 0.1) # get hypomethylated probes
#' Hyper.probe <- get.diff.meth(data, 
#'                             diff.dir="hyper",
#'                             group.col = "definition", 
#'                             sig.dif = 0.1) # get hypomethylated probes
get.diff.meth <- function(data,
                          diff.dir = "hypo",
                          cores = 1,
                          minSubgroupFrac = 0.2,
                          pvalue = 0.01,
                          group.col,
                          min.samples = 5,
                          group1,
                          group2,
                          test = t.test,
                          sig.dif = 0.3,
                          dir.out = "./",
                          save = TRUE){
  
  if(is.null(getMet(data)))
    stop("Cannot identify differential DNA methylation region without DNA methylation data.")
  if(nrow(colData(data))==0){
    stop("Sample information data to do differential analysis.")
  } else if (missing(group.col)){
    stop("Please colData.col should be specified, labeling two group of sample for comparison. See colnames(colData(data)) for possibilities")
  } else if (!group.col %in% colnames(colData(data))){
    stop("Group column not found in phenotypic data and meta-data of the object. See values with colData(data)")
  } else if (missing(group1) | missing(group2)) {
    if(length(unique(colData(data)[,group.col])) < 2){
      stop("Group column should have at least 2 distinct group labels for comparison.")
    } else if (length(unique(colData(data)[,group.col])) > 2){
      stop("Please your object must have only two groups. We found more than two and this might impact the next analysis steps.")
    } else {
      # TO be changed
      groups <- colData(data)[,group.col]
      group1 <- unique(groups)[1] 
      group2 <- unique(groups)[2]
      message(paste0("Group 1: ", group1, "\nGroup 2: ", group2))
    }
  } else if(!group1 %in% unique(colData(data)[,group.col])){
    stop(group1," not found in ", group.col)
  } else if(!group2 %in% unique(colData(data)[,group.col])){
    stop(group2," not found in ", group.col)
  }
  counts <- plyr::count(MultiAssayExperiment::colData(data)[,group.col])
  message(paste0("ELMER will search for probes ", diff.dir,"methylated in group ",
                 group1, " (n:",subset(counts,x == group1)$freq,")", 
                 " compared to ", 
                 group2, " (n:",subset(counts,x == group2)$freq,")"))
  message(paste0("Number of probes: ",nrow(getMet(data))))
  
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  Top.m <- ifelse(diff.dir == "hyper",TRUE,FALSE)
  groups.info <- colData(data)[getMetSamples(data),group.col]
  met <- assay(getMet(data))
  probes <- rownames(met)
  out <- alply(.data = met, .margins = 1,
               .fun = function(x) {
                 Stat.diff.meth(percentage = minSubgroupFrac,
                                meth = x,
                                min.samples = min.samples,
                                groups = groups.info,
                                group1 = group1,
                                test = test,
                                group2 = group2,
                                Top.m = Top.m)},
               .progress = "text", .parallel = parallel
  )
  out <- do.call(rbind,out)
  out <- as.data.frame(out,stringsAsFactors = FALSE)
  out$probe <- probes
  diffCol <- paste0(gsub("[[:punct:]]| ", ".", group1),"_Minus_",gsub("[[:punct:]]| ", ".", group2))
  out$adjust.p <- p.adjust(as.numeric(out$PP),method = "BH")
  out <- out[,c("probe","PP","MeanDiff","adjust.p")]
  colnames(out) <- c("probe","pvalue", diffCol, "adjust.p")
  rownames(out) <- out$probe
  
  if(save){
    message("Saving results")
    dir.create(dir.out,showWarnings = FALSE, recursive = TRUE)
    write_csv(out,path=sprintf("%s/getMethdiff.%s.probes.csv",dir.out,diff.dir))
    write_csv(out[out$adjust.p < pvalue & abs(out[,diffCol]) > sig.dif & !is.na(out$adjust.p),],
              path=sprintf("%s/getMethdiff.%s.probes.significant.csv",dir.out,diff.dir))
  }
  
  result <- out[out$adjust.p < pvalue & abs(out[,diffCol]) > sig.dif & !is.na(out$adjust.p),]
  if(nrow(result) == 0 ) {
    message("No relevant probes found")
  } else {
    message(paste0("Number of relevant probes found: ", nrow(result)))
  }  
  return(result)  
}

## TCGA pipe don't specify dir.out
#' get.pair to predict enhancer-gene linkages.
#' @description 
#' get.pair is a function to predict enhancer-gene linkages using associations between 
#' DNA methylation at enhancer CpG sites and expression of 20 nearby genes of the CpG sites
#' (see reference). Two files will be saved if save is true: getPair.XX.all.pairs.statistic.csv
#' and getPair.XX.pairs.significant.csv (see detail).
#' @usage 
#' get.pair(data, 
#'          nearGenes, 
#'          minSubgroupFrac = 0.4, 
#'          permu.size = 10000, 
#'          permu.dir = NULL, 
#'          raw.pvalue = 0.001, 
#'          Pe = 0.001, 
#'          mode = "unsupervised",
#'          diff.dir = NULL,
#'          dir.out = "./",
#'          diffExp = FALSE,
#'          group.col, 
#'          group1, 
#'          group2, 
#'          cores = 1, 
#'          filter.probes = TRUE, 
#'          filter.portion = 0.3,  
#'          filter.percentage = 0.05,
#'          label = NULL, save = TRUE)
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. 
#' See \code{\link{createMAE}} function.
#' @param nearGenes Can be either a list containing output of GetNearGenes 
#' function or path of rda file containing output of GetNearGenes function.
#' @param cores A interger which defines number of core to be used in parallel process.
#'  Default is 1: don't use parallel process.
#' @param minSubgroupFrac A number ranging from 0 to 1, specifying the fraction of 
#' extreme  samples that define group U (unmethylated) and group M (methylated), 
#' which are used to link probes to genes. 
#' The default is 0.4 (the lowest quintile of samples is the U group and the highest quintile samples is the M group) 
#' because we typically want to be able to detect a specific (possibly unknown) molecular subtype among tumor; 
#' these subtypes often make up only a minority of samples, and 20\% was chosen as a lower bound for the purposes of statistical power. 
#' If you are using pre-defined group labels, such as treated replicates vs. untreated replicated, use a value of 1.0 (Supervised mode).
#' @param permu.size A number specify the times of permuation. Default is 10000.
#' @param permu.dir A path where the output of permutation will be. 
#' @param raw.pvalue A number specify the raw p-value cutoff for defining signficant pairs.
#'  Default is 0.001. It will select the significant P value  cutoff before calculating the empirical p-values.
#' @param Pe A number specify the empirical p-value cutoff for defining signficant pairs.
#'  Default is 0.001
#' @param filter.probes Should filter probes by selecting only probes that have at least
#' a certain number of samples below and above a certain cut-off. 
#' See \code{\link{preAssociationProbeFiltering}} function.
#' @param filter.portion A number specify the cut point to define binary methlation level for probe loci. 
#' Default is 0.3. When beta value is above 0.3, the probe is methylated and 
#' vice versa. For one probe, the percentage of methylated and unmethylated samples 
#' should be above filter.percentage value.   
#' Only used if filter.probes is TRUE. See \code{\link{preAssociationProbeFiltering}} function.
#' @param filter.percentage Minimun percentage of samples to be considered in methylated and unmethylated
#' for the filter.portion option. Default 5\%. Only used if filter.probes is TRUE.
#'  See \code{\link{preAssociationProbeFiltering}} function.
#' @param diffExp A logic. Default is FALSE. If TRUE, t test will be applied to 
#'  test whether putative target gene are differentially expressed between two groups.
#' @param group.col A column defining the groups of the sample. You can view the 
#' available columns using: colnames(MultiAssayExperiment::colData(data)).
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param mode A character. Can be "unsupervised" or "supervised". If unsupervised is set
#' the U (unmethylated) and M (methylated) groups will be selected 
#' among all samples based on methylation of each probe.
#' Otherwise U group and M group will set as the samples of group1 or group2 as described below:
#' If diff.dir is "hypo, U will be the group 1 and M the group2.
#' If diff.dir is "hyper" M group will be the group1 and U the group2.
#' @param diff.dir A character can be "hypo" or "hyper", showing differential 
#' methylation dirction in group 1.  It can be "hypo" which means the probes are hypomethylated in group1; 
#' "hyper" which means the probes are hypermethylated in group1; 
#' This argument is used only when mode is supervised nad 
#' it should be the same value from get.diff.meth function.
#' @param dir.out A path specify the directory for outputs. Default is current directory
#' @param label A character labels the outputs.
#' @param save Two files will be saved if save is true: getPair.XX.all.pairs.statistic.csv
#' and getPair.XX.pairs.significant.csv (see detail).
#' @return Statistics for all pairs and significant pairs
#' @export 
#' @author 
#' Lijing Yao (creator: lijingya@usc.edu) 
#' Tiago C Silva (maintainer: tiagochst@usp.br)
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' data <- ELMER:::getdata("elmer.data.example")
#' nearGenes <-GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
#'                          geneAnnot=getExp(data))
#' Hypo.pair <- get.pair(data=data,
#'                        nearGenes=nearGenes,
#'                        permu.size=5,
#'                        group.col = "definition",
#'                        group1 = "Primary solid Tumor", 
#'                        group2 = "Solid Tissue Normal",
#'                        raw.pvalue = 0.2,
#'                        Pe = 0.2,
#'                        dir.out="./",
#'                        label= "hypo")
#'                      
#' Hypo.pair <- get.pair(data = data,
#'                       nearGenes = nearGenes,
#'                       permu.size = 5,
#'                       raw.pvalue = 0.2,
#'                       Pe = 0.2,
#'                       dir.out = "./",
#'                       diffExp = TRUE, 
#'                       group.col = "definition",
#'                       group1 = "Primary solid Tumor", 
#'                       group2 = "Solid Tissue Normal",
#'                       label = "hypo")                    
get.pair <- function(data,
                     nearGenes,
                     minSubgroupFrac = 0.4,
                     permu.size = 10000,
                     permu.dir = NULL, 
                     raw.pvalue = 0.001,
                     Pe = 0.001,
                     mode = "unsupervised",
                     diff.dir = NULL,
                     dir.out = "./",
                     diffExp = FALSE,
                     group.col,
                     group1 = NULL,
                     group2 = NULL,
                     cores = 1,
                     filter.probes = TRUE,
                     filter.portion = 0.3, 
                     filter.percentage = 0.05,
                     label = NULL,
                     save = TRUE){
  
  if(!all(names(nearGenes) %in% rownames(getMet(data))))
    stop("Probes option should be subset of rownames of methylation matrix.")
  if(is.character(nearGenes)){
    nearGenes <- get(load(nearGenes))
  } else if(!is.list(nearGenes)){
    stop("nearGene option must be a list containing output of GetNearGenes function 
         or path of rda file containing output of GetNearGenes function.")
  }
  if(diffExp & missing(group.col)) 
    stop("Please set group.col argument to test whether putative target gene are differentially expressed between two groups.")
  
  if(missing(group.col)) stop("Please set group.col argument")
  if(missing(group1)) stop("Please set group1 argument")
  if(missing(group2)) stop("Please set group2 argument")
  data <- data[,colData(data)[,group.col] %in% c(group1, group2)]
  
  # Supervised groups
  unmethylated <- methylated <- NULL
  if(mode == "supervised"){
    if(is.null(diff.dir)) stop("For supervised mode please set diff.dir argument (same from the get.diff.meth)")
    if(diff.dir == "hypo"){
      message("Using pre-defined groups. U (unmethylated): ",group1,", M (methylated): ", group2)
      unmethylated <-  which(colData(data)[,group.col]  == group1)
      methylated <-  which(colData(data)[,group.col]  == group2)
    } else {
      message("Using pre-defined groups. U (unmethylated): ",group2,", M (methylated): ", group1)
      unmethylated <-  which(colData(data)[,group.col]  == group2)
      methylated <-  which(colData(data)[,group.col]  == group1)
    }
  } else {
    message("Selecting U (unmethylated) and M (methylated) groups. Each groups has ", minSubgroupFrac * 50,"% of samples")
  }
  # Paralellization code
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  
  data <- preAssociationProbeFiltering(data, K = filter.portion, percentage = filter.percentage)
  
  met <- assay(getMet(data))
  # Probes that were removed from the last steps cannot be verified
  nearGenes <- nearGenes[names(nearGenes) %in% rownames(met)] 
  
  if(length(nearGenes) == 0) {
    message("No probes passed the preAssociationProbeFiltering filter")
    return(NULL)
  }
  exp <- assay(getExp(data))
  message("Calculating Pp (probe - gene) for all nearby genes")
  Probe.gene <- adply(.data = names(nearGenes), .margins = 1,
                      .fun = function(x) {
                        Stat.nonpara(Probe = x,
                                     Meths = met[x,], 
                                     NearGenes = nearGenes,
                                     methy = methylated,
                                     unmethy = unmethylated,
                                     Top = minSubgroupFrac/2, # Each group will have half of the samples
                                     Exps = exp)},
                      .progress = "text", .parallel = parallel, .id = NULL
  )

  rownames(Probe.gene) <- paste0(Probe.gene$Probe,".",Probe.gene$GeneID)
  Probe.gene <- Probe.gene[!is.na(Probe.gene$Raw.p),]
  
  if(save) {
    dir.create(dir.out, showWarnings = FALSE)
    file <- sprintf("%s/getPair.%s.all.pairs.statistic.csv",dir.out, label)
    write_csv(Probe.gene,path=file)
    message(paste("File created:", file))
  }
  
  Probe.gene <- Probe.gene[Probe.gene$Raw.p < raw.pvalue,]
  Probe.gene <- Probe.gene[order(Probe.gene$Raw.p),]
  selected <- Probe.gene
  if(nrow(selected) == 0) {
    message(paste("No significant pairs were found for pvalue =", raw.pvalue))
    return(selected)
  }
  
  #   Probe.gene$logRaw.p <- -log10(Probe.gene$Raw.p)
  GeneID <- unique(Probe.gene[,"GeneID"])
  message(paste("Calculating Pr (random probe - gene). Permutating ", permu.size, "probes for",  length(GeneID), "genes"))
  # get permutation
  permu <- get.permu(data,
                     geneID     = GeneID, 
                     percentage = minSubgroupFrac / 2,
                     rm.probes  = names(nearGenes), 
                     methy      = methylated,
                     unmethy    = unmethylated,
                     permu.size = permu.size, 
                     permu.dir  = permu.dir,
                     cores      = cores)
  # Get empirical p-value
  Probe.gene.Pe <- Get.Pvalue.p(Probe.gene,permu)
  
  if(save) write_csv(Probe.gene.Pe, path=sprintf("%s/getPair.%s.pairs.statistic.with.empirical.pvalue.csv",dir.out, label))
  selected <- Probe.gene.Pe[Probe.gene.Pe$Pe < Pe & !is.na(Probe.gene.Pe$Pe),]
  
  # Change distance from gene to nearest TSS
  selected$Distance <- NULL
  selected <- addDistNearestTSS(data, NearGenes = selected)
  
  if(diffExp){
    message("Calculating differential expression between two groups")
    Exp <- assay(getExp(data)[unique(selected$GeneID),])
    groups <- unique(colData(data)[,group.col])
    prefix <- paste(gsub("[[:punct:]]| ", ".", groups),collapse =  ".vs.")
    log.col <- paste0("log2FC_",prefix)
    diff.col <- paste0(prefix,".diff.pvalue")
    idx1 <- colData(data)[,group.col] == groups[1] 
    idx2 <- colData(data)[,group.col] == groups[2] 
    out <- adply(.data = split(Exp,rownames(Exp)), .margins = 1,
                 .fun = function(x) {
                   test <- t.test(x = x[idx1],y = x[idx2])
                   out <- data.frame("log2FC" = test$estimate[1] - test$estimate[2],
                                     "diff.pvalue" = test$p.value)
                 },
                 .progress = "text", .parallel = parallel,.id = "GeneID"
    )
    add <- out[match(selected$GeneID, out$GeneID),c("log2FC","diff.pvalue")]
    colnames(add) <- c(log.col,diff.col)
    selected <- cbind(selected, add)                                                         
  }
  if(save) write_csv(selected,path=sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, label))
  invisible(gc())
  return(selected)
}

### permutation
#permu.size can be all which mean all the usable probes.
#' get.permu to generate permutation results for calculation of empirical P values for
#' each enhancer-gene linkage.
#' @description 
#' get.permu is a function to use the same statistic model to calculate random enhancer-gene 
#' pairs. Based on the permutation value, empirical P value can be calculated for the 
#' real enhancer-gene pair (see reference).
#' @usage 
#' get.permu(data, 
#'           geneID, 
#'           methy = NULL,
#'           unmethy = NULL,
#'           percentage = 0.2, 
#'           rm.probes = NULL, 
#'           permu.size = 10000, 
#'           permu.dir = NULL, 
#'           cores = 1)
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#' @param geneID A vector lists the genes' ID.
#' @param rm.probes A vector lists the probes name.
#' @param cores A interger which defines number of core to be used in parallel process.
#'  Default is 1: don't use parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples of group 1 and group 2
#' groups used to link probes to genes. Default is 0.2.
#' @param permu.size A number specify the times of permuation. Default is 10000.
#' @param permu.dir A path where the output of permuation will be. 
#' @param methy Index of M (methylated) group.
#' @param unmethy Index of U (unmethylated) group.
#' @return Permutations
#' @importFrom plyr alply
#' @importFrom doParallel registerDoParallel
#' @author 
#' Lijing Yao (creator: lijingya@usc.edu) 
#' Tiago C Silva (maintainer: tiagochst@usp.br)
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @note 
#' Permutation is the most time consuming step. It is recommended to use multiple  
#' cores for this step. Default permutation time is 1000 which may need 12 hrs by 4 cores. 
#' However 10,000 permutations is recommended to get high confidence results. But it may cost 2 days.
#' @export 
#' @examples
#' data <- ELMER:::getdata("elmer.data.example")
#' permu <-get.permu(data = data,
#'                   geneID=rownames(getExp(data)),
#'                   rm.probes=c("cg00329272","cg10097755"),
#'                   permu.size=5)
get.permu <- function(data, 
                      geneID, 
                      methy = NULL,
                      unmethy = NULL,
                      percentage = 0.2, 
                      rm.probes = NULL,
                      permu.size = 10000, 
                      permu.dir = NULL,
                      cores = 1){
  
  ## get usable probes
  usable.probes <- names(getMet(data))
  usable.probes <- usable.probes[!usable.probes %in% rm.probes]
  if(length(usable.probes) < permu.size) 
    stop(sprintf("There is no enough usable probes to perform %s time permutation, 
                 set a smaller permu.size.",permu.size))
  if(!is.numeric(permu.size)) permu.size <- length(usable.probes) 
  
  # Desire for reproducible results
  set.seed(200)
  probes.permu <- sample(usable.probes, size = permu.size, replace = FALSE)
  
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  
  
  # We have two cases to consider:
  # 1) Permutation was not done before
  # 2) It was done before
  # 2.a) We have more probes to evaluate
  # 2.b) We have more genes to evaluate
  # 2.c) More genes and more probes
  # 2.d) No more genes or probes
  # For 1) just do for all genes and probes
  # For 2 a-c do it for new probes, then do for new genes for all probes
  # For 2.d just subset
  permu <- NULL
  tmp.probes <- probes.permu
  tmp.genes <- geneID
  missing.genes <- NULL
  # Check if it isCase 2: Permutation already done
  file <- file.path(permu.dir,"permu.rda")
  if (!is.null(permu.dir)) {
    if (file.exists(file)) {
      temp.space <- new.env()
      permu.file <- get(load(file, temp.space), temp.space)
      rm(temp.space)
      # Does the probe really exists ?
      permu.file <- permu.file[,colnames(permu.file) %in% rownames(getMet(data))]
      tmp.probes <- probes.permu[!probes.permu %in%  colnames(permu.file)]
      if(!all(geneID %in% rownames(permu.file))) { 
        tmp.genes <- rownames(permu.file)
        missing.genes <- geneID[!geneID %in% tmp.genes]
      }
    }
  }  
  permu.meth <- assay(getMet(data)[tmp.probes,,drop=FALSE] )
  exp.data <- assay(getExp(data))
  
  # Should Exps=exp.data[geneID,] to improve performance ?
  # in that case For a second run we will need to look if gene is in the matrix and also probe
  if(length(tmp.probes) > 0) {
    exps <- exp.data[tmp.genes,,drop=FALSE]
    permu <- alply(.data = tmp.probes, .margins = 1,
                   .fun = function(x) {
                     Stat.nonpara.permu(
                       Probe = x,
                       Meths = permu.meth[x,],
                       Gene  = tmp.genes,
                       methy = methy,
                       unmethy = unmethy,
                       Top   = percentage,
                       Exps  = exps)},
                   .progress = "text", .parallel = parallel
    )
    
    permu <- sapply(permu,
                    function(x,geneID){ 
                      x <- x[match(geneID,x[,1]),2]
                    },
                    geneID=tmp.genes,simplify=FALSE)
    
    permu <- do.call(cbind,permu)
    rownames(permu) <- tmp.genes
    colnames(permu) <- tmp.probes
  } 
  
  if(!is.null(permu) & length(file) > 0) {
    if(file.exists(file)){
      # Put genes in the same order before rbind it
      permu.file <- permu.file[match(rownames(permu),rownames(permu.file)),,drop=FALSE]
      permu <- cbind(permu, permu.file)
    }
  } else if(is.null(permu) & length(file) > 0) {
    permu <- permu.file
  }
  
  # For the missing genes calculate for all probes
  if(length(missing.genes) > 0) {
    # Get all probes
    permu.meth <- assay(getMet(data)[colnames(permu),] )
    exps <- exp.data[missing.genes,,drop=FALSE]
    permu.genes <- alply(.data = colnames(permu), .margins = 1,
                         .fun = function(x) {
                           Stat.nonpara.permu(
                             Probe = x,
                             Meths = permu.meth[x,],
                             Gene  = missing.genes,
                             Top   = percentage,
                             methy = methy,
                             unmethy = unmethy,
                             Exps  = exps)},
                         .progress = "text", .parallel = parallel
    )
    
    permu.genes <- sapply(permu.genes,
                          function(x,geneID){ 
                            x <- x[match(geneID,x[,1]),2]
                          },
                          geneID=missing.genes,simplify=FALSE)
    
    permu.genes <- do.call(cbind,permu.genes)
    rownames(permu.genes) <- missing.genes
    colnames(permu.genes) <- colnames(permu)
    # Adding new genes
    # Make sure probes are in the same order
    permu.genes <- permu.genes[,match(colnames(permu.genes),colnames(permu))]
    permu <- rbind(permu,permu.genes)
  }
  
  if(!is.null(permu.dir) & !is.null(permu)) {
    dir.create(permu.dir, showWarnings = FALSE, recursive = TRUE)
    save(permu,file = file.path(permu.dir,"permu.rda"), compress = "xz")
  }
  permu <- permu[geneID,probes.permu, drop = FALSE]
  return(permu)
}

#'promoterMeth
#' @title
#' promoterMeth to calculate associations of gene expression with DNA methylation
#' at promoter regions
#' @description 
#' promoterMeth is a function to calculate associations of gene expression with DNA methylation
#' at promoter regions.
#' @usage 
#' promoterMeth(data, sig.pvalue = 0.01, minSubgroupFrac = 0.4, 
#'              upstream = 200,  downstream = 2000, save = TRUE, cores = 1)
#'@param data A Multi Assay Experiment object with DNA methylation and 
#' gene expression Summarized Experiment objects
#'@param sig.pvalue A number specifies significant cutoff for gene silenced by promoter
#' methylation. Default is 0.01. P value is raw P value without adjustment.
#' @param minSubgroupFrac A number ranging from 0 to 1 
#' specifying the percentage of samples used to create the groups U (unmethylated) 
#' and M (methylated) used to link probes to genes. 
#' Default is 0.4 (lowest quintile of all samples will be in the 
#' U group and the highest quintile of all samples in the M group).
#' @param upstream Number of bp upstream of TSS to consider as promoter region
#' @param downstream  Number of bp downstream of TSS to consider as promoter region
#'@param cores Number of cores to be used in paralellization. Default 1 (no paralellization)
#' @param save A logic. If it is true, the result will be saved.  
#' @importFrom GenomicRanges promoters
#' @importFrom utils write.csv
#' @return A data frame contains genes whose expression significantly anti-correlated
#' with promoter methylation.
#' @examples 
#' \dontrun{
#'   data(elmer.data.example.promoter)
#'   Gene.promoter <- promoterMeth(mae.promoter) 
#' }
#' @export
promoterMeth <- function(data,
                         sig.pvalue = 0.01,
                         minSubgroupFrac = 0.4,
                         upstream = 200,
                         downstream = 2000,
                         save = TRUE,
                         cores = 1){
  
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  
  message("Calculating associations of gene expression with DNA methylation at promoter regions")
  TSS_2K <- promoters(rowRanges(getExp(data)), upstream = upstream, downstream = downstream)
  probes <- rowRanges(getMet(data))
  overlap <- findOverlaps(probes, TSS_2K)
  df <- data.frame(Probe=as.character(names(probes)[queryHits(overlap)]), 
                   GeneID=TSS_2K$ensembl_gene_id[subjectHits(overlap)], stringsAsFactors=FALSE)
  if(nrow(df)==0){
    out <- data.frame(GeneID=c(), Symbol=c(), Raw.p= c())
  } else {
    df <- unique(df)
    ProbeInTSS <- split(df$Probe,df$GeneID)
    message("Calculating average DNA methylation for probes near the same TSS")
    ## calculate average methylation of promoter
    met <- assay(getMet(data))
    Gene.promoter <- lapply(ProbeInTSS, 
                            function(x, METH){
                              meth <- METH[x,]
                              if(length(x)>1){
                                meth <- colMeans(meth,na.rm=TRUE)
                              }  
                              return(meth)
                            },   
                            METH=met)
    
    
    Gene.promoter <- do.call(rbind, Gene.promoter)
    ## make fake NearGene 
    Fake <- data.frame(Symbol = values(getExp(data))[values(getExp(data))$ensembl_gene_id %in% rownames(Gene.promoter),"external_gene_name"],
                       GeneID = values(getExp(data))[values(getExp(data))$ensembl_gene_id %in% rownames(Gene.promoter),"ensembl_gene_id"],
                       Distance = 1,
                       Side = 1, stringsAsFactors=FALSE)
    Fake <- split(Fake, Fake$GeneID)
    exps <- assay(getExp(data))
    
    message("Calculating Pp (probe - gene) for all nearby genes")
    out <- adply(.data = rownames(Gene.promoter)[1:3], .margins = 1,
                 .fun = function(x) {
                   Stat.nonpara(Probe = x,
                                Meths = Gene.promoter[x,], 
                                NearGenes = Fake,
                                Top = minSubgroupFrac/2,
                                Exps = exps)},
                 .progress = "text", .parallel = parallel, .id = NULL
    )
    
    out <- out[,c("GeneID","Symbol","Raw.p")]
    if(save) write.csv(out, 
                       file = "Genes_all_anticorrelated_promoter_methylation.csv",
                       row.names = FALSE)
    out <- out[out$Raw.p < sig.pvalue & !is.na(out$Raw.p),]
  }
  if(nrow(out) == 0) message("No assossiation was found")
  if(save) write.csv(out, 
                     file = "Genes_significant_anticorrelated_promoter_methylation.csv",
                     row.names = FALSE)
  return(out)
}
#' get.enriched.motif to identify the overrepresented motifs in a set of probes (HM450K) regions.
#' @description 
#' get.enriched.motif is a function make use of Probes.motif data from \pkg{ELMER.data}  
#' package to calculate the motif enrichment Odds Ratio and  95\% confidence interval for
#' a given set of probes using fisher test function, after performing the Fisher’s exact test, 
#' the results for all transcription factors are corrected for multiple testing with the Benjamini-Hochberg procedure. 
#' If save is TURE, two output files will be saved: 
#' getMotif.XX.enriched.motifs.rda and getMotif.XX.motif.enrichment.csv (see detail).
#' @usage 
#' get.enriched.motif(data, probes.motif, probes, min.motif.quality = "DS",
#'                    background.probes, lower.OR = 1.1, min.incidence = 10, 
#'                    dir.out = "./", label = NULL, save = TRUE, plot.title=NULL)
#' @param data A multi Assay Experiment from  \code{\link{createMAE}} function.
#' If set and probes.motif/background probes are missing this will be used to get 
#' this other two arguments correctly. This argument is not require, you can set probes.motif and 
#' the backaground.probes manually.
#' @param probes.motif A matrix contains motifs occurrence within probes regions. Probes.motif in 
#' \pkg{ELMER.data} will be used if probes.motif is missing (detail see Probes.motif.hg19.450K in ELMER.data).
#' @param probes A vector lists the name of probes to define the set of probes in which motif enrichment
#' OR and confidence interval will be calculated.
#' @param background.probes A vector lists name of probes which are considered as 
#' background for motif.enrichment  calculation (see detail).
#' @param lower.OR A number specifies the smallest lower boundary of 95\% confidence interval for Odds Ratio.
#' The motif with higher lower boudnary of 95\% confidence interval for Odds Ratio than the number 
#' are the significantly enriched motifs (detail see reference).
#' @param min.incidence A non-negative integer specifies the minimum incidence of motif in the given probes set. 
#' 10 is default.
#' @param pvalue FDR P-value cut off (default 0.05)
#' @param min.motif.quality Minimum motif quality score to consider. 
#' Possible valules: A, B, C , D, AS (A and S), BS (A, B and S), CS (A, B , C and S), DS (all - default) 
#' Description: Each PWM has a quality rating from A to D where 
#' A represents motifs with the highest confidence, and D motifs only weakly describe the pattern with a 
#' limited applications for quantitative analyses. 
#' Special S quality marks the single-box motifs (secondary motif). 
#' Source: http://hocomoco.autosome.ru/help#description_quality_score
#' More information: \url{http://nar.oxfordjournals.org/content/44/D1/D116.full#sec-8}
#' @param dir.out A path. Specifies the directory for outputs. Default is current directory
#' @param label A character. Labels the outputs such as "hypo", "hyper"
#' @param save If save is TURE, two files will be saved: getMotif.XX.enriched.motifs.rda and 
#' getMotif.XX.motif.enrichment.csv (see detail).
#' @param plot.title Plot title. Default: no title.
#' @return A list contains enriched motifs with the probes regions harboring the motif.
#' @export 
#' @details 
#'   background.probes:
#'   For enhancer study, it is better to use probes within distal enhancer probes as
#'   background.probes. For promoter study, it is better to use probes within promoter
#'   regions as background.probes. Because enhancer and promoter have different CG content
#'   and harbors different clusters of TFs motif.
#'   
#'   save:
#'   if save is TRUE, two files will be save on the disk. The first file is 
#'   getMotif.XX.motif.enrichment.csv (XX depends on option label). This file reports 
#'   the Odds Ratio and 95\% confidence interval for these Odds Ratios which pass the 
#'   signficant cutoff (lower.OR and min.incidence). The second file is 
#'   getMotif.XX.enriched.motifs.rda (XX depends on option lable). This file contains
#'   a list R object with enriched motifs as name and probes containing the enriched 
#'   motif as contents. This object will be used in \code{\link{get.TFs}} function.
#'   if save is FALSE, the function will return a R object which is the same with second file.
#'@return A list (R object) with enriched motifs as name and probes containing the enriched 
#'   motif as contents. And hypo.motif.enrichment.pdf plot will be generated.
#' @author 
#' Lijing Yao (creator: lijingya@usc.edu) 
#' @importFrom magrittr divide_by multiply_by %>% add
#' @importFrom plyr alply
#' @importFrom utils data read.csv
#' @importFrom S4Vectors metadata
#' @importFrom dplyr filter
#' @importFrom Matrix colMeans colSums
#' @import ELMER.data
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' probes <- c("cg00329272","cg10097755","cg08928189", "cg17153775","cg21156590",
#' "cg19749688","cg12590404","cg24517858","cg00329272","cg09010107",
#' "cg15386853", "cg10097755", "cg09247779","cg09181054","cg19371916")
#'   data <- tryCatch(ELMER:::getdata("elmer.data.example"), error = function(e) {
#'   message(e)
#'   data(elmer.data.example, envir = environment())
#'   })
#' bg <- rownames(getMet(data))
#' data(Probes.motif.hg38.450K,package = "ELMER.data")
#' enriched.motif <- get.enriched.motif(probes.motif = Probes.motif.hg38.450K,
#'                                      probes = probes,
#'                                      background.probes = bg,
#'                                      pvalue = 1, 
#'                                      min.incidence = 2, 
#'                                      label = "hypo")
#' # If the MAE is set, the background and the probes.motif will 
#' # be automatically set                                     
#' enriched.motif <- get.enriched.motif(data = data,
#'                                      min.motif.quality = "DS",
#'                                      probes=probes,
#'                                      pvalue = 1,
#'                                      min.incidence=2, 
#'                                      label="hypo")
get.enriched.motif <- function(data,
                               probes.motif, 
                               probes,
                               min.motif.quality = "DS",
                               background.probes,
                               pvalue = 0.05,
                               lower.OR = 1.1,
                               min.incidence = 10, 
                               dir.out="./",
                               label = NULL,
                               save = TRUE,
                               plot.title = NULL){
  if(missing(probes.motif)){
    if(missing(data)) stop("Please set probes.motif argument. See ELMER data")
    file <- paste0("Probes.motif.",metadata(data)$genome,".",metadata(data)$met.platform)
    message("Loading object: ",file)
    if(file == "Probes.motif.hg38.450K") probes.motif <- getdata("Probes.motif.hg38.450K")
    if(file == "Probes.motif.hg19.450K") probes.motif <- getdata("Probes.motif.hg19.450K")
    if(file == "Probes.motif.hg38.EPIC") probes.motif <- getdata("Probes.motif.hg38.EPIC")
    if(file == "Probes.motif.hg19.EPIC") probes.motif <- getdata("Probes.motif.hg19.EPIC")
  }  
  all.probes.TF <- probes.motif
  ## here need to be add motif search part.
  if(missing(probes)) stop("probes option should be specified.")
  if(missing(background.probes)){
    if(!missing(data)) {
      background.probes <- as.character(names(getMet(data)))
    } else if(file.exists(sprintf("%s/probeInfo_feature_distal.rda",dir.out))){
      background.probes <- get(load(sprintf("%s/probeInfo_feature_distal.rda",dir.out)))
      background.probes <- as.character(names(background.probes))
    } else {
      message("backaground.probes argument is missing. We will use all probes as background, ", 
              "but for enhancer study, it is better to use probes within distal enhancer probes as background.probes.")
      background.probes <- rownames(all.probes.TF)
    }
  }
  probes <- unique(probes) # A probe should not be considered more than one time
  if(length(probes) < min.incidence) stop("Number of unique prober is smaller than the min.incidence required")
  background.probes <- background.probes[background.probes %in% rownames(all.probes.TF)]
  bg.probes.TF <- all.probes.TF[background.probes,]
  bg.Probes.TF.percent <- Matrix::colMeans(bg.probes.TF) # This is equal to: c/(c+d)
  
  ## load probes for enriched motif ----------------------------------------------
  probes.TF <- all.probes.TF[rownames(all.probes.TF) %in% probes,]
  probes.TF.num <- Matrix::colSums(probes.TF, na.rm=TRUE)
  
  # Odds ratio
  #      p/(1-p)     p * (1-P)   where p = a/(a + b) probes with motif
  # OR =--------- = -----------   where P = c/(c + d) bg probes with motif (entire enhancer probe set)
  #      P/(1-P)     P * (1-p)
  p <- Matrix::colMeans(probes.TF)
  P <- bg.Probes.TF.percent
  sub.enrich.TF <- multiply_by(p,(1-P)) %>%  divide_by(P)  %>%  divide_by(1-p)  
  # Extreme cases: p = 1(likely)/0 (likely) or P = 1 (unlikely) / 0 (likely) 
  # case 1 p:1,P=1 OR = 1/0/1/0 = Inf
  # case 2 p:0,P=0 OR = 0/1/0/1 = NaN
  # case 3 p:1,P=0 OR = 1/0/0/1 = Inf
  # case 4 p:0,P=1 OR = 0/1/1/0 = NaN
  # Cases with NaN p = 0, so we will set OR to 0
  sub.enrich.TF[is.nan(sub.enrich.TF)] <- 0 
  # SD = sqrt(1/a + 1/b + 1/c + 1/d)
  # a is the number of probes within the selected probe set that contain one or more motif occurrences; 
  # b is the number of probes within the selected probe set that do not contain a motif occurrence; 
  # c and d are the same counts within the entire enhancer probe set (background)
  # lower boundary of 95% conf idence interval = exp (ln OR - SD)
  a <- Matrix::colSums(probes.TF)
  b <- nrow(probes.TF) - Matrix::colSums(probes.TF)
  c <- Matrix::colSums(bg.probes.TF )
  d <- nrow(bg.probes.TF) - Matrix::colSums(bg.probes.TF)
  fisher <- plyr::adply(seq_along(1:length(a)),.margins = 1, .fun = function(i)  { 
    x <- fisher.test(matrix(c(a[i],b[i],c[i],d[i]),nrow = 2,ncol = 2))
    ret <- data.frame(x$conf.int[1],x$conf.int[2],x$estimate,x$p.value)
    colnames(ret) <- c("lowerOR","upperOR","OR","p.value")
    ret
  },.id = NULL,.progress = "text")
  rownames(fisher) <- names(a)
  Summary <- data.frame(motif  =  names(a),
                        NumOfProbes = probes.TF.num,
                        PercentageOfProbes = probes.TF.num/length(probes),
                        fisher,
                        FDR = p.adjust(fisher$p.value,method = "BH"),
                        stringsAsFactors = FALSE)
  hocomoco <- getHocomocoTable()
  family.class <- hocomoco[,c("Model",grep("family",colnames(hocomoco),value = T))]
  Summary <- merge(Summary,family.class, by.x = "motif",by.y = "Model")
  Summary <- Summary[order(Summary$lowerOR, decreasing = TRUE),]
  if(save) write_csv(Summary, 
                     path = sprintf("%s/getMotif.%s.motif.enrichment.csv",
                                    dir.out,label))
  
  ## enriched motif and probes
  en.motifs <- as.character(Summary[Summary$lowerOR > lower.OR & Summary$NumOfProbes > min.incidence & Summary$FDR <= pvalue,"motif"])
  
  # Subset by quality
  print.header("Filtering motifs based on quality", "subsection")
  message("Number of enriched motifs with quality:")
  message("-----------")
  for(q in c("A","B","C","D","S"))  message(paste0(" => ",q,": ", length(grep(paste0("\\.",q),en.motifs))))
  message("-----------")
  
  en.motifs <- grep(paste0("\\.[A-",toupper(min.motif.quality),"]"), en.motifs, value = T)
  message("Considering only motifs with quality from A up to ", min.motif.quality,": ",length(en.motifs)," motifs are enriched.")
  enriched.motif <- alply(en.motifs, 
                          function(x, probes.TF) {
                            rownames(probes.TF[probes.TF[,x]==1,x,drop=FALSE])
                          },
                          probes.TF=probes.TF,.margins = 1, .dims = FALSE)
  attributes(enriched.motif) <- NULL
  names(enriched.motif) <- en.motifs
  
  if(save) save(enriched.motif, file = sprintf("%s/getMotif.%s.enriched.motifs.rda",dir.out,label))
  
  ## make plot 
  suppressWarnings({
    P <- motif.enrichment.plot(motif.enrichment = filter(Summary,grepl(paste0("\\.[A-",toupper(min.motif.quality),"]"), Summary$motif)), 
                          significant = list(NumOfProbes = 10, lowerOR = 1.1, OR = 1.3), 
                          dir.out = dir.out,
                          label=paste0(label,".quality.A-",toupper(min.motif.quality)),
                          save=TRUE)
    P <- motif.enrichment.plot(motif.enrichment = filter(Summary,grepl(paste0("\\.[A-",toupper(min.motif.quality),"]"), Summary$motif)), 
                          significant = list(OR = 1.3), 
                          dir.out = dir.out,
                          summary = TRUE,
                          label = paste0(label,".quality.A-",toupper(min.motif.quality),"_with_summary"),
                          title = plot.title,
                          save = TRUE)
  })
  ## add information to siginificant pairs
  if(file.exists(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, label))){
    sig.Pairs <- read.csv(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, label), 
                          stringsAsFactors=FALSE)
    if(all(sig.Pairs$Probe %in% rownames(probes.TF))){
      motif.Info <- sapply(sig.Pairs$Probe,
                           function(x, probes.TF,en.motifs){
                             TFs <- names(probes.TF[x,probes.TF[x,]==1])
                             non.en.motif <- paste(setdiff(TFs,en.motifs),collapse = ";")
                             en.motif <- paste(intersect(TFs,en.motifs), collapse = ";")
                             out <- data.frame(non_enriched_motifs=non.en.motif, 
                                               enriched_motifs=en.motif, 
                                               stringsAsFactors = FALSE)
                             return(out)
                           },
                           probes.TF=probes.TF, en.motifs=en.motifs,simplify=FALSE)
      motif.Info <- do.call(rbind,motif.Info)
      sig.Pairs <- cbind(sig.Pairs, motif.Info)
      write_csv(sig.Pairs, 
                path=sprintf("%s/getPair.%s.pairs.significant.withmotif.csv",dir.out, label))
    }
  }
  
  return(enriched.motif)
}

#' get.TFs to identify regulatory TFs.
#' @description 
#' get.TFs is a function to identify regulatory TFs based on motif analysis and association analysis 
#' between the probes containing a particular motif and expression of all known TFs. If save is true, 
#' two files will be saved: getTF.XX.significant.TFs.with.motif.summary.csv and getTF.hypo.TFs.with.motif.pvalue.rda (see detail).
#' @usage 
#'   get.TFs(data, 
#'           enriched.motif, 
#'           TFs, 
#'           group.col,
#'           group1,
#'           group2,
#'           mode = "unsupervised",
#'          diff.dir = NULL,
#'           motif.relevant.TFs, 
#'           minSubgroupFrac = 0.4,
#'           dir.out = "./",
#'           label = NULL, 
#'           cores = 1,
#'           save = TRUE)
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE}} function.
#' @param enriched.motif A list containing output of get.enriched.motif function or a path of XX.rda file containing output of get.enriched.motif function.
#' @param TFs A data.frame containing TF GeneID and Symbol or a path of XX.csv file containing TF GeneID and Symbol.
#' If missing, human.TF list will be used (human.TF data in ELMER.data). 
#' For detail information, refer the reference paper.
#' @param motif.relevant.TFs A list containing motif as names and relavent TFs as contents
#'  for each list element or a path of XX.rda file containing a list as above. 
#' If missing, motif.relavent.TFs will be used (motif.relavent.TFs data in ELMER.data). 
#' For detail information, refer the reference paper.
#' @param minSubgroupFrac A number ranging from 0 to 1 
#' specifying the percentage of samples used to create the groups U (unmethylated) 
#' and M (methylated) used to link probes to TF expression. 
#' Default is 0.4 (lowest quintile of all samples will be in the 
#' U group and the highest quintile of all samples in the M group).
#' @param cores A interger which defines the number of cores to be used in parallel process. Default is 1: no parallel process.
#' @param dir.out A path specifies the directory for outputs of get.pair function. Default is current directory
#' @param label A character labels the outputs.
#' @param save A logic. If save is ture, two files will be saved: getTF.XX.significant.TFs.with.motif.summary.csv and 
#' getTF.hypo.TFs.with.motif.pvalue.rda (see detail). If save is false, a data frame contains the same content with the first file.
#' @param group.col A column defining the groups of the sample. You can view the 
#' available columns using: colnames(MultiAssayExperiment::colData(data)).
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param mode A character. Can be "unsupervised" or "supervised". If unsupervised is set
#' the U (unmethylated) and M (methylated) groups will be selected 
#' among all samples based on methylation of each probe.
#' Otherwise U group and M group will set as the samples of group1 or group2 as described below:
#' If diff.dir is "hypo, U will be the group 1 and M the group2.
#' If diff.dir is "hyper" M group will be the group1 and U the group2.
#' @param diff.dir A character can be "hypo" or "hyper", showing differential 
#' methylation dirction in group 1.  It can be "hypo" which means the probes are hypomethylated in group1; 
#' "hyper" which means the probes are hypermethylated in group1; 
#' This argument is used only when mode is supervised nad 
#' it should be the same value from get.diff.meth function.
#' @export 
#' @details 
#' save: If save is ture, two files will be saved. The first file is getTF.XX.significant.TFs.with.motif.summary.csv (XX depends on option lable). 
#' This file contain the regulatory TF significantly associate with average DNA methylation at particular motif sites. 
#' The second file is getTF.hypo.TFs.with.motif.pvalue.rda (XX depends on option label). 
#' This file contains a matrix storing the statistic results for significant associations between TFs (row) and average DNA methylation at motifs (column). 
#' If save is false, a data frame which contains the same content with the first file will be reported.
#' @importFrom plyr ldply  adply
#' @importFrom doParallel registerDoParallel
#' @importFrom stats na.omit
#' @importFrom parallel detectCores
#' @return 
#'  Potential responsible TFs will be reported in a dataframe with 4 columns:
#'  \itemize{
#'    \item{motif: the names of motif.}
#'    \item{top.potential.TF.family: the highest ranking upstream TFs which are known recognized the motif. First item in potential.TFs.family}
#'    \item{top.potential.TF.subfamily: the highest ranking upstream TFs which are known recognized the motif. First item in potential.TFs.subfamily}
#'    \item{potential.TFs.family: TFs which are within top 5\% list and are known recognized the motif  (considering family classification).}
#'    \item{potential.TFs.subfamily: TFs which are within top 5\% list and are known recognized the motif (considering subfamily classification).}
#'    \item{top_5percent: all TFs which are within top 5\% list.}
#'  }
#' @author 
#' Lijing Yao (creator: lijingya@usc.edu) 
#' Tiago C Silva (maintainer: tiagochst@usp.br)
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' data <- tryCatch(
#'   ELMER:::getdata("elmer.data.example"), 
#'   error = function(e) {
#'     message(e)
#'      data(elmer.data.example, envir = environment())
#'   })
#' enriched.motif <- list("P53_HUMAN.H11MO.1.A"= c("cg00329272", "cg10097755", "cg08928189",
#'                                  "cg17153775", "cg21156590", "cg19749688", "cg12590404",
#'                                  "cg24517858", "cg00329272", "cg09010107", "cg15386853",
#'                                  "cg10097755", "cg09247779", "cg09181054"))
#' TF <- get.TFs(data, 
#'               enriched.motif, 
#'               group.col = "definition",
#'               group1 = "Primary solid Tumor", 
#'               group2 = "Solid Tissue Normal",
#'               TFs = data.frame(
#'                      external_gene_name=c("TP53","TP63","TP73"),
#'                      ensembl_gene_id= c("ENSG00000141510",
#'                                         "ENSG00000073282",
#'                                         "ENSG00000078900"),
#'                      stringsAsFactors = FALSE),
#'              label="hypo")
#' # This case will use Uniprot dabase to get list of Trasncription factors              
#' TF <- get.TFs(data, 
#'               group.col = "definition",
#'               group1 = "Primary solid Tumor", 
#'               group2 = "Solid Tissue Normal",
#'               enriched.motif, 
#'               label="hypo")              
get.TFs <- function(data, 
                    enriched.motif, 
                    TFs, 
                    group.col,
                    group1,
                    group2,
                    mode = "unsupervised",
                    diff.dir = NULL,
                    motif.relevant.TFs,
                    minSubgroupFrac = 0.4,
                    dir.out = "./",
                    label = NULL,
                    cores = 1,
                    save = TRUE){
  if(missing(data)){
    stop("data argument is empty")
  }
  if(missing(enriched.motif)){
    stop("enriched.motif is empty.")
  }else if(is.character(enriched.motif)){
    enriched.motif <- get(load(enriched.motif)) # The data is in the one and only variable
  } else if(!is.list(enriched.motif)){
    stop("enriched.motif option should be a list object.")
  }
  if(length(enriched.motif) == 0) {
    message("No enriched motifs were found in the last step")
    return(NULL)
  }
  
  if(missing(group.col)) stop("Please set group.col argument")
  if(missing(group1)) stop("Please set group1 argument")
  if(missing(group2)) stop("Please set group2 argument")
  data <- data[,colData(data)[,group.col] %in% c(group1, group2)]
  
  # Supervised groups
  unmethylated <- methylated <- NULL
  if(mode == "supervised"){
    if(is.null(diff.dir)) stop("For supervised mode please set diff.dir argument (same from the get.diff.meth)")
    if(diff.dir == "hypo"){
      message("Using pre-defined groups. U (unmethylated): ",group1,", M (methylated): ", group2)
      unmethylated <-  which(colData(data)[,group.col]  == group1)
      methylated <-  which(colData(data)[,group.col]  == group2)
    } else {
      message("Using pre-defined groups. U (unmethylated): ",group2,", M (methylated): ", group1)
      unmethylated <-  which(colData(data)[,group.col]  == group2)
      methylated <-  which(colData(data)[,group.col]  == group1)
    }
  } else {
    message("Selecting U (unmethylated) and M (methylated) groups. Each groups has ", minSubgroupFrac * 50,"% of samples")
  }
  
  
  if(missing(TFs)){
    # Here we will make some assumptions:
    # TFs has a column Symbol
    # data has the field external_gene_name which should be created by 
    # createMultAssayExperiment function
    # external_gene_name is the default for hg38 in biomart
    # external_gene_id is the default for hg19 in biomart
    TFs <- getTF()
  } else if(is.character(TFs)){
    TFs <- read.csv(TFs, stringsAsFactors=FALSE)
  }
  
  if(missing(motif.relevant.TFs)){
    message("Accessing TF families from TFClass database to indentify known potential TF")
    TF.family <-  createMotifRelevantTfs()
    TF.subfamily <-  createMotifRelevantTfs("subfamily")
  } else if(is.character(motif.relevant.TFs)){
    TF.family <- get(load(motif.relevant.TFs)) # The data is in the one and only variable
    TF.subfamily <-  TF.family
  }
  
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  
  # This will calculate the average methylation at all motif-adjacent probes 
  message("Calculating the average methylation at all motif-adjacent probes ")
  
  motif.meth <- ldply(enriched.motif, 
                      function(x,meth){
                        if(length(x)<2) { 
                          return(meth[x,])
                        } else {
                          return(colMeans(meth[x,],na.rm = TRUE))
                        }}, meth = assay(getMet(data))[unique(unlist(enriched.motif)),,drop = FALSE],
                      .progress = "text", .parallel = parallel, .id = "rownames"
  )
  rownames(motif.meth) <- motif.meth$rownames
  motif.meth$rownames <- NULL
  
  # motif.meth matrix 
  # - rows: average methylation at all motif-adjacent probes (rownames will be the motif)
  # - cols: each patient
  
  # rownames are ensemble gene id
  TFs <- TFs[TFs$ensembl_gene_id %in% rownames(getExp(data)),]
  gene <- TFs$ensembl_gene_id
  gene.name <- TFs$external_gene_name # For plotting purposes 
  
  # Definition:
  # M group: 20% of samples with the highest average methylation at all motif-adjacent probes
  # U group: 20% of samples with the lowest 
  
  # The Mann-Whitney U test was used to test 
  # the null hypothesis that overall gene expression in group M was greater or equal 
  # than that in group U.
  message("Performing Mann-Whitney U test")
  
  # For each motif (x) split the Meths object into U and M and evaluate the expression
  # of all TF Exps (obj)
  exps <- assay(getExp(data))[gene,]
  TF.meth.cor <- alply(.data = names(enriched.motif), .margins = 1,
                       .fun = function(x) {
                         Stat.nonpara.permu( 
                           Probe = x,
                           Meths = motif.meth[x,],
                           Gene  = gene,
                           unmethy = unmethylated,
                           methy = methylated,
                           Top   = minSubgroupFrac/2,
                           Exps  = exps)},
                       .progress = "text", .parallel = parallel
  )
  TF.meth.cor <- lapply(TF.meth.cor, function(x){return(x$Raw.p)})
  TF.meth.cor <- do.call(cbind,TF.meth.cor)
  ## check row and col names
  rownames(TF.meth.cor) <- gene.name
  colnames(TF.meth.cor) <- names(enriched.motif)
  TF.meth.cor <- na.omit(TF.meth.cor)
  
  # TF.meth.cor matrix with raw p-value (Pr)
  # - rows: TFs
  # - cols: motifs
  # lower Raw p-values means that TF expression in M group is lower than in 
  # U group. That means, the Unmethylated with more TF expression
  # have a higher correlation.
  
  message("Finding potential TF and known potential TF")
  
  # For each motif evaluate TF
  cor.summary <- adply(colnames(TF.meth.cor), 
                       function(x, TF.meth.cor, motif.relavent.TFs.family,motif.relavent.TFs.subfamily){ 
                         cor <- rownames(TF.meth.cor)[sort(TF.meth.cor[,x],index.return = TRUE)$ix]
                         top <- cor[1:floor(0.05 * nrow(TF.meth.cor))]
                         if (any(top %in% motif.relavent.TFs.family[[x]])) {
                           potential.TF.family <- top[top %in% motif.relavent.TFs.family[[x]]]
                         } else {
                           potential.TF.family <- NA
                         }
                         if (any(top %in% motif.relavent.TFs.subfamily[[x]])) {
                           potential.TF.subfamily <- top[top %in% motif.relavent.TFs.subfamily[[x]]]
                         } else {
                           potential.TF.subfamily <- NA
                         }
                         out <- data.frame("motif" = x,
                                           "top.potential.TF.family" = ifelse(!is.na(potential.TF.family[1]),potential.TF.family[1],NA),
                                           "top.potential.TF.subfamily" = ifelse(!is.na(potential.TF.subfamily[1]),potential.TF.subfamily[1],NA),
                                           "potential.TF.family" = ifelse(!any(sapply(potential.TF.family,is.na)),paste(potential.TF.family, collapse = ";"),NA),
                                           "potential.TF.subfamily" = ifelse(!any(sapply(potential.TF.subfamily,is.na)),paste(potential.TF.subfamily, collapse = ";"),NA),
                                           "top_5percent_TFs" = paste(top,collapse = ";"),
                                           stringsAsFactors = FALSE)
                         return(out)
                       },                                         
                       TF.meth.cor=TF.meth.cor, motif.relavent.TFs.family=TF.family,motif.relavent.TFs.subfamily=TF.subfamily, 
                       .progress = "text", .parallel = parallel,.margins = 1, .id = NULL)
  rownames(cor.summary) <- cor.summary$motif
  
  if(save){
    save(TF.meth.cor, 
         file=sprintf("%s/getTF.%s.TFs.with.motif.pvalue.rda",dir.out=dir.out, label=label))
    write_csv(cor.summary, 
              path=sprintf("%s/getTF.%s.significant.TFs.with.motif.summary.csv",
                           dir.out=dir.out, label=label))
  } 
  
  print.header("Creating plots", "subsection")
  message("TF rank plot highlighting TF in the same family (folder: ", sprintf("%s/TFrankPlot",dir.out),")")
  dir.create(sprintf("%s/TFrankPlot",dir.out), showWarnings = FALSE, recursive = TRUE)
  TF.rank.plot(motif.pvalue = TF.meth.cor, 
               motif        = colnames(TF.meth.cor), 
               dir.out      = sprintf("%s/TFrankPlot",dir.out), 
               save         = TRUE)
  return(cor.summary)
}