##get.distal.en
#' get.feature.probe
#' @importFrom GenomicRanges promoters 
#' @importFrom minfi getAnnotation
#' @description This function selects the probes on HM450K that either overlap 
#' distal biofeatures or TSS promoter. 
#' @param promoter A logical.If TRUE, function will ouput the promoter probes.
#' If FALSE, function will ouput the distal probes overlaping with features. The 
#' default is FALSE.
#' @param feature A GRange object containing biofeature coordinate such as 
#' enhancer coordinates. Default is comprehensive genomic enhancer regions from REMC and FANTOM5.
#' feature option is only usable when promoter option is FALSE.
#' @param TSS A GRange object containing the transcription start site. Default is UCSC gene TSS.
#' @param TSS.range A list specify how to define promoter regions. 
#' Default is upstream =2000bp and downstream=2000bp.
#' @param rm.chr A vector of chromosome need to be remove from probes such as chrX chrY or chrM
#' @return A GRange object containing probes that satisfy selecting critiria.
#' @export 
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
get.feature.probe <- function(feature,TSS,TSS.range=list(upstream=2000,downstream=2000),
                              promoter=FALSE,rm.chr=NULL){
  probe <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, what="Locations")
  probe <- GRanges(seqnames=probe$chr,
                   ranges=IRanges(probe$pos,
                                  width=1,
                                  names=rownames(probe)),
                   strand=probe$strand,
                   name=rownames(probe))
  if(!is.null(rm.chr)) probe <- probe[!as.character(seqnames(probe)) %in% rm.chr]
  if(!promoter){
    if(missing(TSS)){
      newenv <- new.env()
      data("Combined.TSS",package = "ELMER.data",envir=newenv)
      TSS <- get(ls(newenv)[1],envir=newenv)   
    }
    TSS <- suppressWarnings(promoters(TSS,upstream = TSS.range[["upstream"]], 
                                      downstream = TSS.range[["downstream"]]))
    probe <- probe[setdiff(1:length(probe),unique(queryHits(findOverlaps(probe,TSS))))]
    if(missing(feature)){
      newenv <- new.env()
      data("Union.enhancer",package = "ELMER.data",envir=newenv)
      feature <- get(ls(newenv)[1],envir=newenv)   
      probe <- probe[unique(queryHits(findOverlaps(probe,feature)))]  
    }else if(is(feature,"GRange")){             
      probe <- probe[unique(queryHits(findOverlaps(probe,feature)))]
    }else{
      stop("feature is not GRange object.")
    }
  }else{
    if(missing(TSS)){
      TSS <- txs()
    }
    TSS <- suppressWarnings(promoters(TSS,upstream = TSS.range[["upstream"]], 
                                      downstream = TSS.range[["downstream"]]))
    probe <- probe[unique(queryHits(findOverlaps(probe,TSS)))]
  }
  return(probe)
}

## get differential methylated probes-------------------------
## TCGA pipe don't specify dir.out
#' get.diff.meth to identify hypo/hyper-methylated CpG sites on HM450K between control and experimental groups such as normal verus tumor samples.
#' @description 
#' get.diff.meth applys one-way t-test to identify the CpG sites that are significantly 
#' hypo/hyper-methyalated using proportional samples (defined by percentage option) from control 
#' and experimental groups. The P values will be adjusted by Benjamini-Hochberg method. 
#' Option pvalue and sig.dif will be the criteria (cutoff) for selecting significant differentially methylated CpG sites.
#'  If save is TURE, two getMethdiff.XX.csv files will be generated (see detail).
#' @param mee A MEE.data object containing at least meth and probeInfo.
#' @param diff.dir A character can be "hypo" or "hyper", showing dirction DNA methylation changes. If it is "hypo", 
#' get.diff.meth function will identify all significantly hypomethylated CpG sites; 
#' If "hyper", get.diff.meth function will identify all significantly hypoermethylated CpG sites
#' @param cores A interger which defines the number of cores to be used in parallel process. Default is 1: no parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples from control and experimental groups that are used to identify the differential methylation. Default is 0.2.
#' @param pvalue A number specifies the significant P value (adjusted P value by BH) cutoff for selecting significant hypo/hyper-methylated probes. Default is 0.01
#' @param sig.dif A number specifies the smallest DNA methylation difference as a cutoff for selecting significant hypo/hyper-methylated probes. Default is 0.3.
#' @param dir.out A path specify the directory for outputs. Default is is current directory.
#' @param test Statistical test to be used. Options: t.test, wilcox.test
#' @param save A logic. When TRUE, two getMethdiff.XX.csv files will be generated (see detail)
#' @return Statistics for all probes and significant hypo or hyper-methylated probes.
#' @export 
#' @importFrom plyr adply
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' library(MultiAssayExperiment)
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' exp.sample <- DataFrame(mee@sample[,c("exp.ID","TN")])
#' rownames(exp.sample) <- exp.sample$exp.ID; exp.sample$exp.ID <- NULL
#' met.sample <- DataFrame(mee@sample[,c("meth.ID","TN")])
#' rownames(met.sample) <- met.sample$meth.ID; met.sample$meth.ID <- NULL
#' sample <- rbind( exp.sample,met.sample)
#' mee.mae <- createMultiAssayExperiment(exp = mee@exp, met = mee@meth,pData = sample, TCGA =  TRUE)
#' Hypo.probe <- get.diff.meth(mee.mae, diff.dir="hypo",group.col = "TN") # get hypomethylated probes
get.diff.meth <- function(data,
                          diff.dir="hypo",
                          cores=1,
                          percentage=0.2,
                          pvalue=0.01,
                          group.col, 
                          group1,
                          group2,
                          test = wilcox.test,
                          sig.dif=0.3,
                          dir.out="./",
                          save=TRUE){
  if(is.null(getMet(data)))
    stop("Cannot identify differential DNA methylation region without DNA methylation data.")
  if(nrow(pData(data))==0){
    stop("Sample information data to do differential analysis.")
  }else if(missing(group.col)){
    stop("Please pData.col should be specified, labeling two group of sample for comparison. See colnames(pData(data)) for possibilities")
  } else if(missing(group1) | missing(group2)) {
    if(length(unique(pData(data)[,group.col]))<2){
      stop("Group column should have at least 2 distinct group labels for comparison.")
    } else {
      # TO be changed
      groups <- pData(data)[,group.col]
      group1 <- unique(groups)[1] 
      group2 <- unique(groups)[2]
      message(paste0("Group 1: ", group1, "\nGroup 2: ", group2))
    }
  }
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  Top.m <- ifelse(diff.dir == "hyper",TRUE,FALSE)
  out <- adply(.data = rownames(getMet(data)), .margins = 1,
               .fun = function(x) {
                 Stat.diff.meth( probe = x,
                                 percentage = percentage,
                                 meth=assay(getMet(data)),
                                 groups = pData(data)[,group.col],
                                 group1 = group1,
                                 test = test,
                                 group2 = group2,
                                 Top.m=Top.m)},
               .progress = "text", .parallel = parallel
  )
  out[,1] <- NULL
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
  
  result <- out[out$adjust.p < pvalue & abs(out[,diffCol])>sig.dif,]
  if(nrow(result) == 0 ) message("No relevant probes found")
  
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
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMultiAssayExperiment function}}.
#' @param nearGenes Can be either a list containing output of GetNearGenes 
#' function or path of rda file containing output of GetNearGenes function.
#' @param cores A interger which defines number of core to be used in parallel process.
#'  Default is 1: don't use parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of 
#' samples used to link probes to genes. Default is 0.2.
#' @param permu.size A number specify the times of permuation. Default is 1000.
#' @param permu.dir A path where the output of permuation will be. 
#' @param Pe A number specify the empircal pvalue cutoff for defining signficant pairs.
#'  Default is 0.001
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
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' gene.info <- TCGAbiolinks:::get.GRCh.bioMart()
#' gene.info <- gene.info[match(gsub("ID","",rownames(mee@exp)),gene.info$entrezgene),]
#' exp <- mee@exp
#' rownames(exp) <- gene.info$ensembl_gene_id
#' exp <- exp[!is.na(rownames(exp)),]
#' data <- createMultiAssayExperiment(exp = exp, met = mee@meth, TCGA = T, genome = "hg19" )
#' nearGenes <-GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
#'                          geneAnnot=getExp(data))
#' Hypo.pair <-get.pair(data=data,
#'                      nearGenes=nearGenes,
#'                      permu.size=5,
#'                      Pe = 0.2,
#'                      dir.out="./",
#'                      label= "hypo")
get.pair <- function(data,
                     probes, # necessary ?
                     nearGenes,
                     percentage=0.2,
                     permu.size=10000,
                     permu.dir=NULL, 
                     Pe=0.001,
                     dir.out="./",
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
  Probe.gene <- adply(.data = names(nearGenes), .margins = 1,
                      .fun = function(x) {
                        Stat.nonpara(Probe = x,
                                     Meths= assay(getMet(data)[x,]), 
                                     NearGenes=nearGenes,
                                     K=portion,
                                     Top=percentage,
                                     Exps=assay(getExp(data)))},
                      .progress = "text", .parallel = parallel
  )
  Probe.gene[,1] <- NULL
  rownames(Probe.gene) <- paste0(Probe.gene$Probe,".",Probe.gene$GeneID)
  Probe.gene <- Probe.gene[!is.na(Probe.gene$Raw.p),]
  
  #   Probe.gene$logRaw.p <- -log10(Probe.gene$Raw.p)
  GeneID <- unique(Probe.gene[!is.na(Probe.gene$Raw.p),"GeneID"])
  message("Permutation\n")
  # get permutation
  permu <- get.permu(data,
                     geneID=GeneID, 
                     percentage=percentage, 
                     rm.probes=names(nearGenes), 
                     permu.size=permu.size, 
                     portion = portion,
                     permu.dir=permu.dir,
                     cores=cores)
  #get empirical p-value
  message("Calculate empirical P value.\n")
  Probe.gene.Pe <- Get.Pvalue.p(Probe.gene,permu)
  Probe.gene.Pe <- Probe.gene.Pe[order(Probe.gene.Pe$Raw.p),]
  if(save) write.csv(Probe.gene.Pe, 
                     file=sprintf("%s/getPair.%s.all.pairs.statistic.csv",dir.out, label),
                     row.names=FALSE)
  selected <- Probe.gene.Pe[Probe.gene.Pe$Pe < Pe & !is.na(Probe.gene.Pe$Pe),]
  
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
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMultiAssayExperiment function}}.
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
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' gene.info <- TCGAbiolinks:::get.GRCh.bioMart()
#' gene.info <- gene.info[match(gsub("ID","",rownames(mee@exp)),gene.info$entrezgene),]
#' exp <- mee@exp
#' rownames(exp) <- gene.info$ensembl_gene_id
#' exp <- exp[!is.na(rownames(exp)),]
#' data <- createMultiAssayExperiment(exp = exp, met = mee@meth, TCGA = T, genome = "hg19" )
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
  
  # Why a const seed?
  set.seed(200)
  
  ## get usable probes
  binary.m <- rowMeans(Binary(assay(getMet(data)),portion),na.rm = TRUE)
  usable.probes <- names(binary.m[binary.m < 0.95 & binary.m > 0.05 & !is.na(binary.m)])
  usable.probes <- usable.probes[!usable.probes %in% rm.probes]
  if(length(usable.probes) < permu.size) 
    stop(sprintf("There is no enough usable probes to perform %s time permutation, 
                 set a smaller permu.size.",permu.size))
  if(!is.numeric(permu.size)) permu.size <- length(usable.probes) 
  probes.permu <- sample(usable.probes, size = permu.size, replace = FALSE)
  
  if(is.null(permu.dir)){
    permu.meth <- assay(getMet(data)[probes.permu,] )
    
    parallel <- FALSE
    if (cores > 1){
      if (cores > detectCores()) cores <- detectCores()
      registerDoParallel(cores)
      parallel = TRUE
    }
    exp.data <- assay(getExp(data))
    permu <- alply(.data = probes.permu, .margins = 1,
                   .fun = function(x) {
                     Stat.nonpara.permu(
                       Probe = x,
                       Meths=permu.meth[x,],
                       Gene=geneID,
                       Top=percentage,
                       Exps=exp.data[geneID,])},
                   .progress = "text", .parallel = parallel
    )

    permu <- sapply(permu,
                    function(x,geneID){ 
                      x <- x[match(geneID,x[,1]),2]},
                    geneID=geneID,simplify=FALSE)
  } else {
    ## if file already there don't need to calculate.
    if(!file.exists(permu.dir)){
      dir.create(permu.dir,recursive = TRUE)
    }
    if(!all(probes.permu %in% dir(permu.dir))){
      tmp.probes <- probes.permu[!probes.permu %in% dir(permu.dir)]
      permu.meth <-  assay(getMet(data)[tmp.probes,])
      parallel <- FALSE
      if (cores > 1){
        if (cores > detectCores()) cores <- detectCores()
        registerDoParallel(cores)
        parallel = TRUE
      }
      exp.data <- assay(getExp(data))
      permu <- alply(.data = tmp.probes, .margins = 1,
                     .fun = function(x) {
                       Stat.nonpara.permu(
                         Probe = x,
                         Meths=permu.meth[x,],
                         Gene=geneID, # How to do it  
                         Top=percentage,
                         Exps=exp.data,
                         permu.dir=permu.dir)},
                     .progress = "text", .parallel = parallel
      )
      
    }
    permu.p <- paste0(permu.dir,"/",probes.permu)
    permu <- sapply(permu.p,
                    function(x,geneID){ 
                      tmp <- read.table(x,stringsAsFactors=FALSE)
                      tmp <- tmp[match(geneID,tmp[,1]),2]},
                    geneID=geneID,simplify=FALSE)
  }
  permu <- do.call(cbind,permu)
  rownames(permu) <- geneID
  colnames(permu) <- probes.permu
  return(permu)
}

#'promoterMeth
#'@param mee A MEE.data object must contains meth, exp, probeInfo, geneInfo four component.
#'@param sig.pvalue A number specify significant cutoff for gene silenced by promoter
#'methylation. Default is 0.01
#' @param percentage A number ranges from 0 to 1 specifying the percentage of 
#' samples used to link probes to genes. Default is 0.2.
#' @importFrom GenomicRanges promoters
#' @return A data frame contains genes whose expression significantly anti-correlated 
#' with promoter methylation.
#' @export
promoterMeth <- function(data,
                         sig.pvalue=0.01,
                         percentage=0.2,
                         save=TRUE){
  TSS_2K <- promoters(rowRanges(getExp(data)), upstream = 100, downstream = 700)
  probes <- rowRanges(getMet(data))
  overlap <- findOverlaps(probes, TSS_2K)
  df <- data.frame(Probe=as.character(names(probes)[queryHits(overlap)]), 
                   GeneID=TSS_2K$ensembl_gene_id[subjectHits(overlap)], stringsAsFactors=FALSE)
  if(nrow(df)==0){
    out <- data.frame(GeneID=c(), Symbol=c(), Raw.p= c())
  }else{
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
                       Distance= 1,
                       Side = 1, stringsAsFactors=FALSE)
    Fake <- split(Fake, Fake$GeneID)
    out <- lapply(rownames(Gene.promoter),Stat.nonpara, NearGenes=Fake,K=0.3,Top=0.2,
                  Meths=Gene.promoter,Exps=assay(getExp(data)))
    out <- do.call(rbind, out)[,c("GeneID","Symbol","Raw.p")]
    out <- out[out$Raw.p < sig.pvalue & !is.na(out$Raw.p),]
  }
  if(save) write.csv(out, file="Genes_significant_anticorrelated_promoter_methylation.csv",
                     row.names=FALSE)
  return(out)
}

#' get.enriched.motif
#' @param probes.motif A matrix contains motifs occurrence within probes regions. 
#' @param probes A vector lists the probes' names in which motif enrichment will be calculated.
#' @param background.probes A vector list of probes' names which are considered as 
#' background for motif.enrichment calculation.
#' @param lower.OR A number specify the lower boundary of Odds ratio which defines 
#' the significant enriched motif. 1.1 is default.
#' @param min.incidence A non-negative integer specify the minimum incidence of 
#' motif in the given probes set. 10 is default.
#' @param dir.out A path specify the directory for outputs. Default is current directory
#' @param label A character labels the outputs.
#' @return A list contains enriched motifs with the probes regions harboring the motif.
#' @export 
#' @examples
#' probes <- c("cg00329272","cg10097755","cg08928189", "cg17153775","cg21156590",
#' "cg19749688","cg12590404","cg24517858","cg00329272","cg09010107",
#' "cg15386853", "cg10097755", "cg09247779","cg09181054","cg19371916")
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' bg <- rownames(getMeth(mee))
#' enriched.motif <- get.enriched.motif(probes=probes,
#'                                      background.probes = bg,
#'                                      min.incidence=2, label="hypo")
get.enriched.motif <- function(probes.motif, 
                               probes, 
                               background.probes,
                               lower.OR=1.1,
                               min.incidence=10, 
                               dir.out="./",
                               label=NULL,
                               save=TRUE){
  if(missing(probes.motif)){
    newenv <- new.env()
    data("Probes.motif",package = "ELMER.data",envir=newenv)
    all.probes.TF <- get(ls(newenv)[1],envir=newenv) 
    # The data is in the one and only variable
  }else{
    all.probes.TF <- probes.motif
  }
  ## here need to be add motif search part.
  if(missing(probes)) stop("probes option should be specified.")
  if(missing(background.probes)){
    if(file.exists(sprintf("%s/probeInfo_feature_distal.rda",dir.out))){
      background.probes <- get(load(sprintf("%s/probeInfo_feature_distal.rda",dir.out)))
      background.probes <- as.character(background.probes$name)
    }else{
      background.probes <- rownames(all.probes.TF)
    }
  }
  bg.probes.TF <- all.probes.TF[background.probes,]
  bg.Probes.TF.percent <- colMeans(bg.probes.TF)
  ## load probes for enriched motif ----------------------------------------------
  probes.TF <- all.probes.TF[probes,]
  probes.TF.num <- colSums(probes.TF, na.rm=TRUE)
  sub.enrich.TF <- colMeans(probes.TF)*(1-bg.Probes.TF.percent)/bg.Probes.TF.percent/(1-colMeans(probes.TF))
  SE <- sqrt(1/colSums(probes.TF) + 1/(nrow(probes.TF)-colSums(probes.TF)) +
               1/colSums(bg.probes.TF)+ 1/(nrow(bg.probes.TF)-colSums(bg.probes.TF)))
  sub.enrich.TF.lower <- exp(log(sub.enrich.TF)-1.96*SE)
  sub.enrich.TF.upper <- exp(log(sub.enrich.TF)+1.96*SE)
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
  en.motifs <- names(sub.enrich.TF.lower[sub.enrich.TF.lower > lower.OR &
                                           !sub.enrich.TF.lower %in% "Inf" & 
                                           probes.TF.num > min.incidence])
  message(sprintf("%s motifs are enriched.",length(en.motifs)))
  enriched.motif <- sapply(en.motifs, 
                           function(x, probes.TF) {
                             names(probes.TF[probes.TF[,x]==1,x])
                           },
                           probes.TF=probes.TF)
  if(save) save(enriched.motif, file= sprintf("%s/getMotif.%s.enriched.motifs.rda",dir.out,label))
  
  ## make plot----
  motif.enrichment.plot(motif.enrichment=Summary, 
                        significant=list(OR=1.3), dir.out =dir.out,label=label, save=TRUE)
  
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
#' between the probes containing a particular motif and expression of all known TFs. If save is ture, 
#' two files will be saved: getTF.XX.significant.TFs.with.motif.summary.csv and getTF.hypo.TFs.with.motif.pvalue.rda (see detail).
#' @usage get.TFs(mee, enriched.motif, TFs, motif.relavent.TFs, percentage = 0.2, dir.out = "./", label = NULL, cores = NULL,save=TRUE)
#' @param data #' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. See \code{\link{createMultiAssayExperiment function}}.
#' @param enriched.motif A list containing output of get.enriched.motif function or a path of XX.rda file containing output of get.enriched.motif function.
#' @param TFs A data.frame containing TF GeneID and Symbol or a path of XX.csv file containing TF GeneID and Symbol.
#' If missing, human.TF list will be used (human.TF data in ELMER.data). 
#' For detail information, refer the reference paper.
#' @param motif.relavent.TFs A list containing motif as names and relavent TFs as contents
#'  for each list element or a path of XX.rda file containing a list as above. 
#' If missing, motif.relavent.TFs will be used (motif.relavent.TFs data in ELMER.data). 
#' For detail information, refer the reference paper.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples of control and experimental groups used to link probes to genes. Default is 0.2.
#' @param cores A interger which defines the number of cores to be used in parallel process. Default is 1: no parallel process.
#' @param dir.out A path specifies the directory for outputs of get.pair function. Default is current directory
#' @param label A character labels the outputs.
#' @param save A logic. If save is ture, two files will be saved: getTF.XX.significant.TFs.with.motif.summary.csv and 
#' getTF.hypo.TFs.with.motif.pvalue.rda (see detail). If save is false, a data frame contains the same content with the first file.
#' @return Potential responsible TFs will be reported.
#' @export 
#' @details 
#' save: If save is ture, two files will be saved. The first file is getTF.XX.significant.TFs.with.motif.summary.csv (XX depends on option lable). 
#' This file contain the regulatory TF significantly associate with average DNA methylation at particular motif sites. 
#' The second file is getTF.hypo.TFs.with.motif.pvalue.rda (XX depends on option label). 
#' This file contains a matrix storing the statistic results for significant associations between TFs (row) and average DNA methylation at motifs (column). 
#' If save is false, a data frame which contains the same content with the first file will be reported.
#' @importFrom pbapply pbsapply
#' @importFrom plyr alply
#' @importFrom doParallel registerDoParallel
#' @author 
#' Lijing Yao (creator: lijingya@usc.edu) 
#' Tiago C Silva (maintainer: tiagochst@usp.br)
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' @examples
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' gene.info <- TCGAbiolinks:::get.GRCh.bioMart()
#' gene.info <- gene.info[match(gsub("ID","",rownames(mee@exp)),gene.info$entrezgene),]
#' exp <- mee@exp
#' rownames(exp) <- gene.info$ensembl_gene_id
#' exp <- exp[!is.na(rownames(exp)),]
#' data <- createMultiAssayExperiment(exp = exp, met = mee@meth, TCGA = T, genome = "hg19" )
#' enriched.motif <- list("TP53"= c("cg00329272", "cg10097755", "cg08928189",
#'                                  "cg17153775", "cg21156590", "cg19749688", "cg12590404",
#'                                  "cg24517858", "cg00329272", "cg09010107", "cg15386853",
#'                                  "cg10097755", "cg09247779", "cg09181054"))
#' TF <- get.TFs(data, enriched.motif, 
#'               TFs=data.frame(GeneID=c("7157","8626","7161"),
#'               Symbol=c("TP53","TP63","TP73"),
#'                ENSEMBL= c("ENSG00000141510","ENSG00000073282","ENSG00000078900"),
#'               stringsAsFactors = FALSE),
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
    enriched.motif <- get(enriched.motif) # The data is in the one and only variable
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
    Tfs <- getTF()
    
  } else if(is.character(TFs)){
    TFs <- read.csv(TFs, stringsAsFactors=FALSE)
  }
  
  if(missing(motif.relavent.TFs)){
    newenv <- new.env()
    data("motif.relavent.TFs",package = "ELMER.data",envir=newenv)
    motif.relavent.TFs <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  } else if(is.character(motif.relavent.TFs)){
    motif.relavent.TFs <- get(load(motif.relavent.TFs)) # The data is in the one and only variable
  }
  
  motif.meth <- lapply(enriched.motif, 
                       function(x,meth){
                         if(length(x)<2) { 
                           return(meth[x,])
                         } else {
                           return(colMeans(meth[x,],na.rm = TRUE))
                         }}, meth = assay(getMet(data))[unique(unlist(enriched.motif)),])
  
  motif.meth <- do.call(rbind, motif.meth)
  
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  if(all(grepl("ENSG",rownames(getExp(data))))) {
    gene <-  TFs$ENSEMBL
  } else {
    gene <- TFs$Symbol
  }
  TF.meth.cor <- alply(.data = names(enriched.motif), .margins = 1,
                       .fun = function(x) {
                         Stat.nonpara.permu( 
                           Probe = x,
                           Meths=motif.meth,
                           Gene=TFs$ENSEMBL,
                           Top=percentage,
                           Exps=assay(getExp(data)))},
                       .progress = "text", .parallel = parallel
  )
  TF.meth.cor <- lapply(TF.meth.cor, function(x){return(x$Raw.p)})
  TF.meth.cor <- do.call(cbind,TF.meth.cor)
  ## check row and col names
  rownames(TF.meth.cor) <- TFs$Symbol
  colnames(TF.meth.cor) <- names(enriched.motif)
  
  cor.summary <- sapply(colnames(TF.meth.cor), 
                        function(x, TF.meth.cor, motif.relavent.TFs){ 
                          cor <- sort(TF.meth.cor[,x])
                          top <- names(cor[1:floor(0.05*nrow(TF.meth.cor))])
                          potential.TF <- top[top %in% motif.relavent.TFs[[x]]]
                          out <- data.frame("motif" = x,
                                            "top potential TF" = potential.TF[1],
                                            "potential TFs" = paste(potential.TF, collapse = ";"),
                                            "top_5percent" = paste(top,collapse = ";"))
                        },                                         
                        TF.meth.cor=TF.meth.cor, motif.relavent.TFs=motif.relavent.TFs, simplify=FALSE)
  cor.summary <- do.call(rbind, cor.summary)
  if(save){
    save(TF.meth.cor, 
         file=sprintf("%s/getTF.%s.TFs.with.motif.pvalue.rda",dir.out=dir.out, label=label))
    write.csv(cor.summary, 
              file=sprintf("%s/getTF.%s.significant.TFs.with.motif.summary.csv",
                           dir.out=dir.out, label=label), row.names=TRUE)
  } 
  if(!file.exists(sprintf("%s/TFrankPlot",dir.out)))
    dir.create(sprintf("%s/TFrankPlot",dir.out))
  TF.rank.plot(motif.pvalue=TF.meth.cor, motif=colnames(TF.meth.cor), 
               dir.out=sprintf("%s/TFrankPlot",dir.out), save=TRUE)
  return(cor.summary)
}