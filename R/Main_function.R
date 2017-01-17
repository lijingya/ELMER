#' @title get.feature.probe to select probes within promoter regions or distal regions.
#' @description 
#' get.feature.probe is a function to select the probes falling into 
#' distal feature regions or promoter regions.
#' @importFrom GenomicRanges promoters 
#' @importFrom minfi getAnnotation
#' @description This function selects the probes on HM450K that either overlap 
#' distal biofeatures or TSS promoter. 
#' @param promoter A logical.If TRUE, function will ouput the promoter probes.
#' If FALSE, function will ouput the distal probes overlaping with features. The 
#' default is FALSE.
#' @param met.platform DNA methyaltion platform to retrieve data from: EPIC or 450K (default)
#' @param genome Genome  GENCODE
#' 
#' @param feature A GRange object containing biofeature coordinate such as 
#' enhancer coordinates. Default is comprehensive genomic enhancer regions from 
#' REMC and FANTOM5 which is Union.enhancer data in \pkg{ELMER.data}.
#' feature option is only usable when promoter option is FALSE.
#' @param TSS A GRange object contains the transcription start sites. When promoter is FALSE, Union.TSS
#' in \pkg{ELMER.data} will be used for default. When promoter is TRUE, UCSC gene TSS will
#' be used as default (see detail). User can specify their own preference TSS annotation. 
#' @param TSS.range A list specify how to define promoter regions. 
#' Default is upstream =2000bp and downstream=2000bp.
#' @param rm.chr A vector of chromosome need to be remove from probes such as chrX chrY or chrM
#' @return A GRange object containing probes that satisfy selecting critiria.
#' @export 
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
get.feature.probe <- function(feature,
                              TSS,
                              genome = "hg38",
                              met.platform = "450K",
                              TSS.range=list(upstream=2000,downstream=2000),
                              promoter=FALSE,
                              rm.chr=NULL){
  probe <- getInfiniumAnnotation(toupper(met.platform),genome)
  # We will rmeove the rs probes, as they should not be used in the analysis
  probe <- probe[!grepl("rs",names(probe)),]
  probe <- probe[!probe$MASK.mapping,] # remove masked probes
  if(!is.null(rm.chr)) probe <- probe[!as.character(seqnames(probe)) %in% rm.chr]
  if(!promoter){
    if(missing(TSS)){
      # The function getTSS gets the transcription coordinantes from Ensemble (GENCODE)
      TSS.gencode <- getTSS(genome = genome)
      # The function txs gets the transcription coordinantes from USCS
      TSS <- TSS.gencode
      # TSS.uscs <- txs(genome)
      # TSS <- c(TSS.gencode, reduce(TSS.uscs))
    }
    suppressWarnings({
      TSS <- promoters(TSS,upstream = TSS.range[["upstream"]], 
                       downstream = TSS.range[["downstream"]])
      
      probe <- probe[setdiff(1:length(probe),unique(queryHits(findOverlaps(probe,TSS))))]
    })
    
    if(missing(feature)){
      newenv <- new.env()
      if(genome == "hg19") data("Union.enhancer.hg19",package = "ELMER.data", envir = newenv)
      if(genome == "hg38") data("Union.enhancer.hg38",package = "ELMER.data", envir = newenv)
      feature <- get(ls(newenv)[1],envir=newenv)   
      probe <- probe[unique(queryHits(findOverlaps(probe,feature)))]
    } else if(is(feature,"GRange")) {             
      probe <- probe[unique(queryHits(findOverlaps(probe,feature)))]
    } else {
      stop("feature is not GRange object.")
    }
  } else {
    if(missing(TSS)){
      # The function getTSS gets the transcription coordinantes from Ensemble (GENCODE)
      TSS.gencode <- getTSS(genome = genome)
      TSS <- TSS.gencode
      # The function txs gets the transcription coordinantes from USCS
      #TSS.uscs <- txs(genome)
      #TSS <- c(TSS.gencode, reduce(TSS.uscs))
    }
    suppressWarnings({
      promoters <- promoters(TSS,upstream = TSS.range[["upstream"]], 
                             downstream = TSS.range[["downstream"]])
      probe <- probe[unique(queryHits(findOverlaps(probe,promoters)))]
    })
  }
  return(probe)
}


## get differential methylated probes-------------------------
## TCGA pipe don't specify dir.out
#' get.diff.meth to identify hypo/hyper-methylated CpG sites on HM450K between control and experimental 
#' groups such as normal verus tumor samples.
#' @description 
#' get.diff.meth applys one-way t-test to identify the CpG sites that are significantly 
#' hypo/hyper-methyalated using proportional samples (defined by percentage option) from control 
#' and experimental groups. The P values will be adjusted by Benjamini-Hochberg method. 
#' Option pvalue and sig.dif will be the criteria (cutoff) for selecting significant differentially methylated CpG sites.
#'  If save is TURE, two getMethdiff.XX.csv files will be generated (see detail).
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE function}}.
#' @param diff.dir A character can be "hypo" or "hyper", showing dirction DNA methylation changes. If it is "hypo", 
#' get.diff.meth function will identify all significantly hypomethylated CpG sites; 
#' If "hyper", get.diff.meth function will identify all significantly hypoermethylated CpG sites
#' @param cores A interger which defines the number of cores to be used in parallel process. Default is 1: no parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples from control and
#'  experimental groups that are used to identify the differential methylation. Default is 0.2.
#' @param pvalue A number specifies the significant P value (adjusted P value by BH) cutoff 
#' for selecting significant hypo/hyper-methylated probes. Default is 0.01
#' @param sig.dif A number specifies the smallest DNA methylation difference as a cutoff for 
#' selecting significant hypo/hyper-methylated probes. Default is 0.3.
#' @param dir.out A path specify the directory for outputs. Default is is current directory.
#' @param test Statistical test to be used. Options: t.test (DEFAULT), wilcox.test
#' @param save A logic. When TRUE, two getMethdiff.XX.csv files will be generated (see detail)
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
#' @importFrom plyr adply
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' data(elmer.data.example)
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
                          percentage = 0.2,
                          pvalue = 0.01,
                          group.col, 
                          group1,
                          group2,
                          test = t.test,
                          sig.dif = 0.3,
                          dir.out = "./",
                          save = TRUE){
  if(is.null(getMet(data)))
    stop("Cannot identify differential DNA methylation region without DNA methylation data.")
  if(nrow(pData(data))==0){
    stop("Sample information data to do differential analysis.")
  } else if (missing(group.col)){
    stop("Please pData.col should be specified, labeling two group of sample for comparison. See colnames(pData(data)) for possibilities")
  } else if (!group.col %in% colnames(pData(data))){
    stop("Group column not found in phenotypic data and meta-data of the object. See values with pData(data)")
  } else if (missing(group1) | missing(group2)) {
    if(length(unique(pData(data)[,group.col]))<2){
      stop("Group column should have at least 2 distinct group labels for comparison.")
    } else if (length(unique(pData(data)[,group.col])) > 2){
      stop("Please your object must have only two groups. We found more than two and this might impact the next analysis steps.")
    } else {
      # TO be changed
      groups <- pData(data)[,group.col]
      group1 <- unique(groups)[1] 
      group2 <- unique(groups)[2]
      message(paste0("Group 1: ", group1, "\nGroup 2: ", group2))
    }
  }
  
  message(paste0("ELMER will search for probes ", diff.dir,"methylated in group ", group1, " compared to ", group2))
  message(paste0("Number of probes: ",nrow(getMet(data))))
  
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  Top.m <- ifelse(diff.dir == "hyper",TRUE,FALSE)
  out <- alply(.data = rownames(getMet(data)), .margins = 1,
               .fun = function(x) {
                 Stat.diff.meth(probe = x,
                                percentage = percentage,
                                meth = assay(getMet(data))[x,],
                                groups = pData(data)[getMetSamples(data),group.col],
                                group1 = group1,
                                test = test,
                                group2 = group2,
                                Top.m=Top.m)},
               .progress = "text", .parallel = parallel
  )
  out <- do.call(rbind,out)
  out <- as.data.frame(out,stringsAsFactors = FALSE)
  diffCol <- paste0(gsub("[[:punct:]]| ", ".", group1),"_Minus_",gsub("[[:punct:]]| ", ".", group2))
  out$adjust.p <- p.adjust(as.numeric(out[,2]),method="BH")
  colnames(out) <- c("probe","pvalue", diffCol, "adjust.p")
  rownames(out) <- out$probe
  if(save){
    write.csv(out,file=sprintf("%s/getMethdiff.%s.probes.csv",dir.out,diff.dir), row.names=FALSE)
    write.csv(out[out$adjust.p < pvalue & abs(out[,diffCol])>sig.dif,],
              file=sprintf("%s/getMethdiff.%s.probes.significant.csv",dir.out,diff.dir), 
              row.names=FALSE)
  }
  
  result <- out[out$adjust.p < pvalue & abs(out[,diffCol])>sig.dif &  !is.na(out$adjust.p),]
  if(nrow(result) == 0 ) {
    message("No relevant probes found")
  } else {
    message(paste0("Number of relevant probes found:", nrow(result)))
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
#' get.pair(data, probes, nearGenes, percentage = 0.2, permu.size = 10000, permu.dir = NULL,  
#'          Pe = 0.001, dir.out = "./", diffExp = FALSE, cores = 1, portion=0.3,   
#'          label = NULL, save=TRUE)
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE function}}.
#' @param nearGenes Can be either a list containing output of GetNearGenes 
#' function or path of rda file containing output of GetNearGenes function.
#' @param cores A interger which defines number of core to be used in parallel process.
#'  Default is 1: don't use parallel process.
#' @param percentage A number ranges from 0 to 0.5 specifying the percentage of 
#' samples used to link probes to genes. Default is 0.2.
#' @param permu.size A number specify the times of permuation. Default is 10000.
#' @param permu.dir A path where the output of permutation will be. 
#' @param calculate.empirical.p A logic. It TRUE (default) execute the permutation step, if FALSE uses the adjusted raw p-value to select the significant pairs.
#' @param pvalue A number specify the pvalue cutoff for defining signficant pairs.
#'  Default is 0.001. If calculate.empirical.p is TRUE it will be the empirical pvalue cutoff.
#'  Otherwise  it will select the significant P value (adjusted P value by BH) cutoff.
#'  @param portion A number specify the cut point for methylated and unmethylated.
#'  Default is 0.3.
#'  @param diffExp A logic. Default is FALSE. If TRUE, t test will be applied to 
#'  test whether putative target gene are differentially expressed between two groups.
#' @param dir.out A path specify the directory for outputs. Default is current directory
#' @param label A character labels the outputs.
#' @return Statistics for all pairs and significant pairs
#' @export 
#' @author 
#' Lijing Yao (creator: lijingya@usc.edu) 
#' Tiago C Silva (maintainer: tiagochst@usp.br)
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' data(elmer.data.example)
#' nearGenes <-GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
#'                          geneAnnot=getExp(data))
#' Hypo.pair <- get.pair(data=data,
#'                      nearGenes=nearGenes,
#'                      permu.size=5,
#'                      pvalue = 0.2,
#'                      dir.out="./",
#'                      label= "hypo")
get.pair <- function(data,
                     nearGenes,
                     percentage=0.2,
                     permu.size=10000,
                     permu.dir=NULL, 
                     pvalue = 0.001,
                     dir.out="./",
                     calculate.empirical.p = TRUE,
                     diffExp=FALSE,
                     cores=1,
                     portion = 0.3, 
                     label=NULL,
                     save=TRUE){
  ## check data
  if(!all(names(nearGenes) %in% rownames(getMet(data))))
    stop("Probes option should be subset of rownames of methylation matrix.")
  if(is.character(nearGenes)){
    nearGenes <- get(load(nearGenes))
  } else if(!is.list(nearGenes)){
    stop("nearGene option must be a list containing output of GetNearGenes function 
         or path of rda file containing output of GetNearGenes function.")
  }
  #get raw pvalue
  ##I need to modify that if there is all NA. stop the process.
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  message("Calculating Pp (probe - gene) for all nearby genes")
  Probe.gene <- adply(.data = names(nearGenes), .margins = 1,
                      .fun = function(x) {
                        Stat.nonpara(Probe = x,
                                     Meths= assay(getMet(data)[x,]), 
                                     NearGenes=nearGenes,
                                     K=portion,
                                     Top=percentage,
                                     Exps=assay(getExp(data)))},
                      .progress = "text", .parallel = parallel, .id = NULL
  )
  rownames(Probe.gene) <- paste0(Probe.gene$Probe,".",Probe.gene$GeneID)
  Probe.gene <- Probe.gene[!is.na(Probe.gene$Raw.p),]
  
  if(calculate.empirical.p){
    #   Probe.gene$logRaw.p <- -log10(Probe.gene$Raw.p)
    GeneID <- unique(Probe.gene[!is.na(Probe.gene$Raw.p),"GeneID"])
    message(paste("Calculating Pr (random probe - gene). Permutating ", permu.size, "probes for each nearby gene"))
    # get permutation
    permu <- get.permu(data,
                       geneID     = GeneID, 
                       percentage = percentage, 
                       rm.probes  = names(nearGenes), 
                       permu.size = permu.size, 
                       portion    = portion,
                       permu.dir  = permu.dir,
                       cores      = cores)
    # Get empirical p-value
    Probe.gene.Pe <- Get.Pvalue.p(Probe.gene,permu)
    
    Probe.gene.Pe <- Probe.gene.Pe[order(Probe.gene.Pe$Raw.p),]
    if(save) write.csv(Probe.gene.Pe, 
                       file=sprintf("%s/getPair.%s.all.pairs.statistic.csv",dir.out, label),
                       row.names = FALSE)
    selected <- Probe.gene.Pe[Probe.gene.Pe$Pe < pvalue & !is.na(Probe.gene.Pe$Pe),]
  } else {
    Probe.gene$Raw.p.adjust <- p.adjust(as.numeric(Probe.gene$Raw.p),method="BH")
    Probe.gene <- Probe.gene[order(Probe.gene$Raw.p.adjust),]
    if(save) write.csv(Probe.gene, 
                       file=sprintf("%s/getPair.%s.all.pairs.statistic.csv",dir.out, label),
                       row.names=FALSE)
    selected <- Probe.gene[Probe.gene$Raw.p.adjust < pvalue & !is.na(Probe.gene$Raw.p.adjust),]
  }
  if(diffExp){
    ## calculate differential expression between two groups.
    Exp <- getExp(mee, geneID = unique(selected$GeneID))
    TN <- getSample(mee,cols = "TN")
    out <- lapply(split(Exp,rownames(Exp)),
                  function(x, TN){
                    test <- t.test(x~TN)
                    U.test <- wilcox.test(x~TN)
                    out <- data.frame("log2FC_TvsN" = test$estimate[2]-test$estimate[1],
                                      "TN.diff.pvalue"=test$p.value)
                    return(out)}, TN=TN)
    out <- do.call(rbind, out)
    out$GeneID <- rownames(out)
    add <- out[match(selected$GeneID, out$GeneID),c("log2FC_TvsN","TN.diff.pvalue")]
    selected <- cbind(selected, add)                                                         
  }
  if(save) write.csv(selected, 
                     file=sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, label),
                     row.names=FALSE)
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
#' get.permu(data, geneID, percentage = 0.2, rm.probes = NULL, portion = 0.3, permu.size = 10000, permu.dir = NULL, cores = 1)
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE function}}.
#' @param geneID A vector lists the genes' ID.
#' @param rm.probes A vector lists the probes name.
#' @param cores A interger which defines number of core to be used in parallel process.
#'  Default is 1: don't use parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples of group 1 and group 2
#' groups used to link probes to genes. Default is 0.2.
#' @param permu.size A number specify the times of permuation. Default is 10000.
#' @param permu.dir A path where the output of permuation will be. 
#' @param portion A number specify the cut point to define binary methlation level for probe loci. 
#' Default is 0.3. When beta value is above 0.3, the probe is methylated and 
#' vice versa. For one probe, the percentage of methylated or unmethylated samples 
#' should be above 0.05.
#' Default is 0.3.
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
#' data(elmer.data.example)
#' permu <-get.permu(data = data,
#'                   geneID=rownames(getExp(data)),
#'                   rm.probes=c("cg00329272","cg10097755"),
#'                   permu.size=5)
get.permu <- function(data, 
                      geneID, 
                      percentage=0.2, 
                      rm.probes=NULL,
                      portion=0.3,
                      permu.size=10000, 
                      permu.dir=NULL,
                      cores=1){
  
  ## get usable probes
  binary.m <- rowMeans(assay(getMet(data)) > portion, na.rm = TRUE)
  usable.probes <- names(binary.m[binary.m < 0.95 & binary.m > 0.05 & !is.na(binary.m)])
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
  # Case 2) Permutation already done
  file <- file.path(permu.dir,"permu.rda")
  if (!is.null(permu.dir)) {
    if (file.exists(file)) {
      temp.space <- new.env()
      permu.file <- get(load(file, temp.space), temp.space)
      rm(temp.space)
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
    permu <- alply(.data = tmp.probes, .margins = 1,
                   .fun = function(x) {
                     Stat.nonpara.permu(
                       Probe = x,
                       Meths=permu.meth[x,],
                       Gene=tmp.genes,
                       Top=percentage,
                       Exps=exp.data[tmp.genes,,drop=FALSE])},
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
      permu.file <- permu.file[match(rownames(permu),rownames(permu.file)),]
      permu <- cbind(permu, permu.file)
    }
  }
  
  # For the missing genes calculate for all probes
  if(length(missing.genes) > 0) {
    # Get all probes
    permu.meth <- assay(getMet(data)[colnames(permu),] )
    permu.genes <- alply(.data = colnames(permu), .margins = 1,
                         .fun = function(x) {
                           Stat.nonpara.permu(
                             Probe = x,
                             Meths=permu.meth[x,],
                             Gene=missing.genes,
                             Top=percentage,
                             Exps=exp.data[missing.genes,,drop=FALSE])},
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
    permu <- permu[match(rownames(permu),geneID),
                   match(colnames(permu),probes.permu)]
  }
  
  if(!is.null(permu.dir) & !is.null(permu)) {
    dir.create(permu.dir, showWarnings = FALSE, recursive = TRUE)
    save(permu,file = file.path(permu.dir,"permu.rda"), compress = "xz")
  }
  if(is.null(permu)) {
    permu <- permu.file[geneID,probes.permu]
  }
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
#' promoterMeth(mee, sig.pvalue = 0.01, percentage = 0.2, save = TRUE)
#'@param mee A Multi Assay Experiment object with DNA methylation and 
#' gene expression Summarized Experiment objects
#'@param sig.pvalue A number specifies significant cutoff for gene silenced by promoter
#' methylation. Default is 0.01. P value is raw P value without adjustment.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of
#'  samples of control and experimental groups used to link promoter DNA methylation to genes.
#'  Default is 0.2.
#' @param save A logic. If it is true, the result will be saved.  
#' @importFrom GenomicRanges promoters
#' @return A data frame contains genes whose expression significantly anti-correlated
#' with promoter methylation.
#' @examples 
#' data(elmer.data.example)
#' Gene.promoter <- promoterMeth(data) 
#' @export
promoterMeth <- function(data,
                         sig.pvalue = 0.01,
                         percentage = 0.2,
                         save = TRUE){
  message("Calculating associations of gene expression with DNA methylation at promoter regions")
  TSS_2K <- promoters(rowRanges(getExp(data)), upstream = 100, downstream = 700)
  probes <- rowRanges(getMet(data))
  overlap <- findOverlaps(probes, TSS_2K)
  df <- data.frame(Probe=as.character(names(probes)[queryHits(overlap)]), 
                   GeneID=TSS_2K$ensembl_gene_id[subjectHits(overlap)], stringsAsFactors=FALSE)
  if(nrow(df)==0){
    out <- data.frame(GeneID=c(), Symbol=c(), Raw.p= c())
  } else {
    df <- unique(df)
    ProbeInTSS <- split(df$Probe,df$GeneID)
    
    ## calculate average methylation of promoter
    Gene.promoter <- lapply(ProbeInTSS, 
                            function(x, METH){
                              meth <- METH[x,]
                              if(length(x)>1){
                                meth <- colMeans(meth,na.rm=TRUE)
                              }  
                              return(meth)
                            },   
                            METH=assay(getMet(data)))
    Gene.promoter <- do.call(rbind, Gene.promoter)
    ## make fake NearGene 
    Fake <- data.frame(Symbol = values(getExp(data))[values(getExp(data))$ensembl_gene_id %in% rownames(Gene.promoter),"external_gene_name"],
                       GeneID = rownames(Gene.promoter),
                       Distance = 1,
                       Side = 1, stringsAsFactors=FALSE)
    Fake <- split(Fake, Fake$GeneID)
    out <- lapply(rownames(Gene.promoter),Stat.nonpara, NearGenes=Fake,K=0.3,Top=0.2,
                  Meths=Gene.promoter,Exps=assay(getExp(data)))
    out <- do.call(rbind, out)[,c("GeneID","Symbol","Raw.p")]
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
#' a given set of probes. If save is TURE, two output files will be saved: 
#' getMotif.XX.enriched.motifs.rda and getMotif.XX.motif.enrichment.csv (see detail).
#' @usage 
#' get.enriched.motif(probes.motif, probes, 
#'                    background.probes, lower.OR = 1.1, min.incidence = 10, 
#'                    dir.out = "./", label = NULL, save=TRUE)
#' @param data A multi Assay Experiment from  \code{\link{createMAE}} function.
#' If set and probes.motif/background probes are missing this will be used to get 
#' this other two arguments correctly. This argument is not require, you can set probes.motif and 
#' the backaground.probes manually.
#' @param probes.motif A matrix contains motifs occurrence within probes regions. Probes.motif in 
#' \pkg{ELMER.data} will be used if probes.motif is missing (detail see \code{\link{Probes.motif}}).
#' @param probes A vector lists the name of probes to define the set of probes in which motif enrichment
#' OR and confidence interval will be calculated.
#' @param background.probes A vector lists name of probes which are considered as 
#' background for motif.enrichment  calculation (see detail).
#' @param lower.OR A number specifies the smallest lower boundary of 95\% confidence interval for Odds Ratio.
#' The motif with higher lower boudnary of 95\% confidence interval for Odds Ratio than the number 
#' are the significantly enriched motifs (detail see reference).
#' @param min.incidence A non-negative integer specifies the minimum incidence of motif in the given probes set. 
#' 10 is default.
#' @param min.motif.quality Minimum motif quality score to consider. 
#' Possible valules: A, B (default), C , D, AS (A and S), BS (A, B and S), CS (A, B , C and S), DS (all) 
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
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' probes <- c("cg00329272","cg10097755","cg08928189", "cg17153775","cg21156590",
#' "cg19749688","cg12590404","cg24517858","cg00329272","cg09010107",
#' "cg15386853", "cg10097755", "cg09247779","cg09181054","cg19371916")
#' data(elmer.data.example)
#' bg <- rownames(getMet(data))
#' enriched.motif <- get.enriched.motif(probes.motif=Probes.motif.hg38.450K,
#'                                      probes=probes,
#'                                      background.probes = bg,
#'                                      min.incidence=2, label="hypo")
#'  # If the MAE is set, the background and the probes.motif will be automatically set                                     
#' enriched.motif <- get.enriched.motif(data = data,
#'                                      min.motif.quality = "DS",
#'                                      probes=probes,
#'                                      min.incidence=2, 
#'                                      label="hypo")
get.enriched.motif <- function(data,
                               probes.motif, 
                               probes,
                               min.motif.quality = "C",
                               background.probes,
                               lower.OR = 1.1,
                               min.incidence = 10, 
                               dir.out="./",
                               label=NULL,
                               save=TRUE){
  if(missing(probes.motif)){
    if(missing(data)) stop("Please set probes.motif argument. See ELMER data")
    file <- paste0("Probes.motif.",metadata(data)$genome,".",metadata(data)$met.platform)
    message("Loading object: ",file)
    newenv <- new.env()
    if(file == "Probes.motif.hg38.450K") data("Probes.motif.hg38.450K", package = "ELMER.data",envir=newenv)
    if(file == "Probes.motif.hg19.450K") data("Probes.motif.hg19.450K", package = "ELMER.data",envir=newenv)
    if(file == "Probes.motif.hg38.EPIC") data("Probes.motif.hg38.EPIC", package = "ELMER.data",envir=newenv)
    if(file == "Probes.motif.hg19.EPIC") data("Probes.motif.hg19.EPIC", package = "ELMER.data",envir=newenv)
    probes.motif <- get(ls(newenv)[1],envir=newenv)   
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
  
  # SD = sqrt(1/a + 1/b + 1/c + 1/d)
  # a is the number of probes within the selected probe set that contain one or more motif occurrences; 
  # b is the number of probes within the selected probe set that do not contain a motif occurrence; 
  # c and d are the same counts within the entire enhancer probe set (background)
  # lower boundary of 95% conf idence interval = exp (ln OR − SD)
  a <- Matrix::colSums(probes.TF)
  b <- nrow(probes.TF) - Matrix::colSums(probes.TF)
  c <- Matrix::colSums(bg.probes.TF )
  d <- nrow(bg.probes.TF) - Matrix::colSums(bg.probes.TF)
  SD <- add(1/a,1/b) %>% add(1/c) %>% add(1/d) %>% sqrt
  sub.enrich.TF.lower <- exp(log(sub.enrich.TF) - 1.96 * SD)
  sub.enrich.TF.upper <- exp(log(sub.enrich.TF) + 1.96 * SD)
  
  ## summary
  Summary <- data.frame(motif = colnames(probes.TF), 
                        NumOfProbes = probes.TF.num,
                        OR = sub.enrich.TF, 
                        lowerOR = sub.enrich.TF.lower, 
                        upperOR = sub.enrich.TF.upper)
  
  Summary <- Summary[order(Summary$lowerOR, decreasing = TRUE),]
  if(save) write.csv(Summary, file= sprintf("%s/getMotif.%s.motif.enrichment.csv",
                                            dir.out,label))
  
  ## enriched motif and probes
  en.motifs <- sub.enrich.TF.lower[sub.enrich.TF.lower > lower.OR &
                                     !sub.enrich.TF.lower %in% "Inf" & 
                                     probes.TF.num > min.incidence]
  # Subset by quality
  print.header("Filtering motifs based on quality", "subsection")
  message("Number of enriched motifs with quality:")
  message("-----------")
  for(q in c("A","B","C","D","S")) message(paste0(" => ",q,": ", length(grep(paste0("H10MO.",q),names(en.motifs)))))
  message("-----------")
  
  en.motifs <- names(en.motifs[grep(paste0("H10MO.[A-",toupper(min.motif.quality),"]"),
                                    names(en.motifs), value = T)])
  message("Considering only motifs with quality from A up to ", min.motif.quality,": ",length(en.motifs)," motifs are enriched.")
  enriched.motif <- alply(en.motifs, 
                          function(x, probes.TF) {
                            rownames(probes.TF[probes.TF[,x]==1,x,drop=FALSE])
                          },
                          probes.TF=probes.TF,.margins = 1, .dims = FALSE)
  attributes(enriched.motif) <- NULL
  names(enriched.motif) <- en.motifs
  
  if(save) save(enriched.motif, file= sprintf("%s/getMotif.%s.enriched.motifs.rda",dir.out,label))
  
  ## make plot 
  suppressWarnings({
    motif.enrichment.plot(motif.enrichment = Summary, 
                          significant = list(OR = 1.3), 
                          dir.out = dir.out,
                          label=label, 
                          save=TRUE)
  })
  ## add information to siginificant pairs
  if(file.exists(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, label))){
    sig.Pairs <- read.csv(sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, label), 
                          stringsAsFactors=FALSE)
    if(all(sig.Pairs$Probe %in% rownames(probes.TF))){
      motif.Info <- sapply(sig.Pairs$Probe,
                           function(x, probes.TF,en.motifs)
                           {TFs <- names(probes.TF[x,probes.TF[x,]==1])
                           non.en.motif <- paste(setdiff(TFs,en.motifs),collapse = ";")
                           en.motif <- paste(intersect(TFs,en.motifs), collapse = ";")
                           out <- data.frame(non_enriched_motifs=non.en.motif, 
                                             enriched_motifs=en.motif, stringsAsFactors = FALSE)
                           return(out)},
                           probes.TF=probes.TF, en.motifs=en.motifs,simplify=FALSE)
      
      motif.Info <- do.call(rbind,motif.Info)
      sig.Pairs <- cbind(sig.Pairs, motif.Info)
      write.csv(sig.Pairs, 
                file=sprintf("%s/getPair.%s.pairs.significant.withmotif.csv",dir.out, label),
                row.names=FALSE)
    }
  }
  return(enriched.motif)
}

#' get.TFs to identify regulatory TFs.
#' @description 
#' get.TFs is a function to identify regulatory TFs based on motif analysis and association analysis 
#' between the probes containing a particular motif and expression of all known TFs. If save is true, 
#' two files will be saved: getTF.XX.significant.TFs.with.motif.summary.csv and getTF.hypo.TFs.with.motif.pvalue.rda (see detail).
#' @usage get.TFs(mee, enriched.motif, TFs, motif.relavent.TFs, percentage = 0.2, dir.out = "./", label = NULL, cores = NULL,save=TRUE)
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMAE function}}.
#' @param enriched.motif A list containing output of get.enriched.motif function or a path of XX.rda file containing output of get.enriched.motif function.
#' @param TFs A data.frame containing TF GeneID and Symbol or a path of XX.csv file containing TF GeneID and Symbol.
#' If missing, human.TF list will be used (human.TF data in ELMER.data). 
#' For detail information, refer the reference paper.
#' @param motif.relavent.TFs A list containing motif as names and relavent TFs as contents
#'  for each list element or a path of XX.rda file containing a list as above. 
#' If missing, motif.relavent.TFs will be used (motif.relavent.TFs data in ELMER.data). 
#' For detail information, refer the reference paper.
#' @param percentage A number ranges from 0 to 0.5 specifying the percentage of samples of control and experimental groups used to link probes to genes. Default is 0.2.
#' @param cores A interger which defines the number of cores to be used in parallel process. Default is 1: no parallel process.
#' @param dir.out A path specifies the directory for outputs of get.pair function. Default is current directory
#' @param label A character labels the outputs.
#' @param save A logic. If save is ture, two files will be saved: getTF.XX.significant.TFs.with.motif.summary.csv and 
#' getTF.hypo.TFs.with.motif.pvalue.rda (see detail). If save is false, a data frame contains the same content with the first file.
#' @export 
#' @details 
#' save: If save is ture, two files will be saved. The first file is getTF.XX.significant.TFs.with.motif.summary.csv (XX depends on option lable). 
#' This file contain the regulatory TF significantly associate with average DNA methylation at particular motif sites. 
#' The second file is getTF.hypo.TFs.with.motif.pvalue.rda (XX depends on option label). 
#' This file contains a matrix storing the statistic results for significant associations between TFs (row) and average DNA methylation at motifs (column). 
#' If save is false, a data frame which contains the same content with the first file will be reported.
#' @importFrom pbapply pbsapply
#' @importFrom plyr ldply  adply
#' @importFrom doParallel registerDoParallel
#' @return 
#'  Potential responsible TFs will be reported in a dataframe with 4 columns:
#'  \itemize{
#'    \item{motif: the names of motif.}
#'    \item{top.potential.TF: the highest ranking upstream TFs which are known recognized the motif. First item in potential.TFs.}
#'    \item{potential.TFs: TFs which are within top 5\% list and are known recognized the motif.}
#'    \item{top_5percent: all TFs which are within top 5\% list.}
#'  }
#' @author 
#' Lijing Yao (creator: lijingya@usc.edu) 
#' Tiago C Silva (maintainer: tiagochst@usp.br)
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' data(elmer.data.example)
#' enriched.motif <- list("P53_HUMAN.H10MO.B"= c("cg00329272", "cg10097755", "cg08928189",
#'                                  "cg17153775", "cg21156590", "cg19749688", "cg12590404",
#'                                  "cg24517858", "cg00329272", "cg09010107", "cg15386853",
#'                                  "cg10097755", "cg09247779", "cg09181054"))
#' TF <- get.TFs(data, 
#'               enriched.motif, 
#'               TFs=data.frame(external_gene_name=c("TP53","TP63","TP73"),
#'                              ensembl_gene_id= c("ENSG00000141510","ENSG00000073282","ENSG00000078900"),
#'                              stringsAsFactors = FALSE),
#'              label="hypo")
#' # This case will use Uniprot dabase to get list of Trasncription factors              
#' TF <- get.TFs(data, 
#'               enriched.motif, 
#'               label="hypo")              
get.TFs <- function(data, 
                    enriched.motif, 
                    TFs, 
                    motif.relavent.TFs,
                    percentage=0.2,
                    dir.out="./",
                    label=NULL,
                    cores=1,
                    save=TRUE){
  if(missing(enriched.motif)){
    stop("enriched.motif is empty.")
  }else if(is.character(enriched.motif)){
    enriched.motif <- get(load(enriched.motif)) # The data is in the one and only variable
  }else if(!is.list(enriched.motif)){
    stop("enriched.motif option should be a list object.")
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
  
  if(missing(motif.relavent.TFs)){
    message("Accessing TF families from TFClass database to indentify known potential TF")
    motif.relavent.TFs <- createMotifRelevantTfs()
  } else if(is.character(motif.relavent.TFs)){
    motif.relavent.TFs <- get(load(motif.relavent.TFs)) # The data is in the one and only variable
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
  TF.meth.cor <- alply(.data = names(enriched.motif), .margins = 1,
                       .fun = function(x) {
                         Stat.nonpara.permu( 
                           Probe = x,
                           Meths=motif.meth[x,],
                           Gene=gene,
                           Top=percentage,
                           Exps=assay(getExp(data))[gene,])},
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
                       function(x, TF.meth.cor, motif.relavent.TFs){ 
                         cor <- rownames(TF.meth.cor)[sort(TF.meth.cor[,x],index.return=T)$ix]
                         top <- cor[1:floor(0.05*nrow(TF.meth.cor))]
                         if (any(top %in% motif.relavent.TFs[[x]])) {
                           potential.TF <- top[top %in% motif.relavent.TFs[[x]]]
                         } else {
                           potential.TF <- NA
                         }
                         out <- data.frame("motif" = x,
                                           "top potential TF" = ifelse(!is.na(potential.TF[1]),potential.TF[1],NA),
                                           "potential TFs" = ifelse(!any(sapply(potential.TF,is.na)),paste(potential.TF, collapse = ";"),NA),
                                           "top_5percent" = paste(top,collapse = ";"),
                                           stringsAsFactors = FALSE)
                         return(out)
                       },                                         
                       TF.meth.cor=TF.meth.cor, motif.relavent.TFs=motif.relavent.TFs, 
                       .progress = "text", .parallel = parallel,.margins = 1, .id = NULL)
  rownames(cor.summary) <- cor.summary$motif
  if(save){
    save(TF.meth.cor, 
         file=sprintf("%s/getTF.%s.TFs.with.motif.pvalue.rda",dir.out=dir.out, label=label))
    write.csv(cor.summary, 
              file=sprintf("%s/getTF.%s.significant.TFs.with.motif.summary.csv",
                           dir.out=dir.out, label=label), row.names=TRUE)
  } 
  if(!file.exists(sprintf("%s/TFrankPlot",dir.out)))
    dir.create(sprintf("%s/TFrankPlot",dir.out))
  print.header("Creating plots", "subsection")
  TF.rank.plot(motif.pvalue=TF.meth.cor, 
               motif=colnames(TF.meth.cor), 
               dir.out=sprintf("%s/TFrankPlot",dir.out), 
               save=TRUE)
  return(cor.summary)
}