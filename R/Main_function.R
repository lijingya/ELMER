##get.distal.en
#' get.feature.probe
#' @import IRanges GenomicRanges
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
#' get.diff.meth
#' @param mee A MEE.data object containing at least meth and probeInfo.
#' @param diff.dir A character can be "hypo" or "hyper", showing differential 
#' methylation dirction.  It can be "hypo" which is only selecting hypomethylated probes; 
#' "hyper" which is only selecting hypermethylated probes; 
#' @param cores A interger which defines number of core to be used in parallel process. 
#' Default is NULL: don't use parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of 
#' samples used to identify the differential methylation. Default is 0.2.
#' @param pvalue A number specify the significant Pvalue cutoff for significant 
#' hypo/hyper-methylated probes. Default is current directory.
#' @param sig.dif A number specify the significant methylation difference cutoff 
#' for significant hypo/hyper-methylated probes. Default is 0.3.
#' @param dir.out A path specify the directory for outputs. Default is is current directory.
#' @param save A logic. When TRUE, output file will be saved.
#' @return Statistics for all probes and significant hypo or hyper-methylated probes.
#' @export 
#' @examples
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' Hypo.probe <- get.diff.meth(mee, diff.dir="hypo") # get hypomethylated probes

get.diff.meth <- function(mee,diff.dir="hypo",cores=NULL,percentage=0.2,
                          pvalue=0.01, sig.dif=0.3, dir.out="./",save=TRUE){
  if(nrow(mee@meth)==0) 
    stop("Cannot identify differential DNA methylation region without DNA methylation data.")
  if(nrow(getSample(mee))==0){
    stop("Sample information data to do differential analysis.")
  }else if(is.null(getSample(mee,cols="TN"))){
    stop("\"TN\" should be specified, labeling two group of sample for comparison.")
  }else if(length(table(getSample(mee,cols="TN")))<2){
    stop("\"TN\" should have at 2 distinct group labels for comparison.")
  }
  
  if(requireNamespace("parallel", quietly=TRUE) && requireNamespace("snow", quietly=TRUE)) {
	  if(!is.null(cores)) {
		  if(cores > parallel::detectCores()) cores <- parallel::detectCores()/2
		  cl <- snow::makeCluster(cores,type = "SOCK")
	  }
  }
  if("hyper" %in% diff.dir){
	  if(requireNamespace("parallel", quietly=TRUE)) {
		  if(!is.null(cores)){
			  out <- parallel::parSapplyLB(cl,rownames(mee@meth),Stat.diff.meth,
										   percentage=percentage,meth=mee@meth,
										   TN=getSample(mee,cols="TN"),Top.m=TRUE,simplify =FALSE)
		  } else {
			  out <- sapply(rownames(mee@meth),Stat.diff.meth,
							percentage=percentage,meth=mee@meth,
							TN=getSample(mee,cols="TN"),Top.m=TRUE,simplify =FALSE) 
		  }
	  } else {
		  out <- sapply(rownames(mee@meth),Stat.diff.meth,
						percentage=percentage,meth=mee@meth,
						TN=getSample(mee,cols="TN"),Top.m=TRUE,simplify =FALSE)
	  }
    out <- do.call(rbind,out)
    out <- as.data.frame(out,stringsAsFactors = FALSE)
    out$adjust.p <- p.adjust(as.numeric(out[,2]),method="BH")
    colnames(out) <- c("probe","pvalue","ExperimentMinControl","adjust.p")
    if(save){
      write.csv(out,file=sprintf("%s/getMethdiff.hyper.probes.csv",dir.out), row.names=FALSE)
      write.csv(out[out$adjust.p < pvalue & abs(out$ExperimentMinControl)>sig.dif,],
                file=sprintf("%s/getMethdiff.hyper.probes.significant.csv",dir.out), 
                row.names=FALSE)
    }
    result <- out[out$adjust.p < pvalue & abs(out$ExperimentMinControl)>sig.dif,]
  }
  if("hypo" %in% diff.dir){
	  if(requireNamespace("parallel")) {
		  if(!is.null(cores)){
			out <- parallel::parSapplyLB(cl,rownames(mee@meth),Stat.diff.meth,
										 percentage=percentage,meth=mee@meth,
										 TN=getSample(mee,cols="TN"),Top.m=FALSE,simplify =FALSE)
			parallel::stopCluster(cl)
		  } else {
			out <- sapply(rownames(mee@meth),Stat.diff.meth,percentage=percentage,
						  meth=mee@meth,TN=getSample(mee,cols="TN"),Top.m=FALSE,
						  simplify =FALSE)
		  }
    }else{
      out <- sapply(rownames(mee@meth),Stat.diff.meth,percentage=percentage,
                    meth=mee@meth,TN=getSample(mee,cols="TN"),Top.m=FALSE,
                    simplify =FALSE)
    }
    out <- do.call(rbind,out)
    out <- as.data.frame(out,stringsAsFactors = FALSE)
    out$adjust.p <- p.adjust(as.numeric(out[,2]),method="BH")
    colnames(out) <- c("probe","pvalue","ExperimentMinControl","adjust.p")
    if(save){
      write.csv(out,file=sprintf("%s/getMethdiff.hypo.probes.csv",dir.out), 
                row.names=FALSE)
      write.csv(out[out$adjust.p < pvalue & abs(out$ExperimentMinControl)>sig.dif,],
                file=sprintf("%s/getMethdiff.hypo.probes.significant.csv",dir.out),
                row.names=FALSE)
    }
    result <- out[out$adjust.p < pvalue & abs(out$ExperimentMinControl)>sig.dif,]
  }
  return(result)  
}


### pairing function--------------------
#nearGenes : can be either a list containing output of GetNearGenes function or 
#path of rda file containing output of GetNearGenes function.
## TCGA pipe don't specify dir.out
#' get pair
#' @param mee A MEE.data object containing at least meth, exp, probeInfo, geneInfo.
#' @param probes A vector lists the probes' name that need to be linked to genes.
#' @param nearGenes Can be either a list containing output of GetNearGenes 
#' function or path of rda file containing output of GetNearGenes function.
#' @param cores A interger which defines number of core to be used in parallel process.
#'  Default is NULL: don't use parallel process.
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
#' @examples
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' nearGenes <-GetNearGenes(TRange=getProbeInfo(mee,probe=c("cg00329272","cg10097755")),
#'                          geneAnnot=getGeneInfo(mee))
#'                          Hypo.pair <-get.pair(mee=mee,probes=c("cg00329272","cg10097755"),
#'                                               nearGenes=nearGenes,permu.size=5,Pe = 0.2,
#'                                               dir.out="./",
#'                                               label= "hypo")
get.pair <- function(mee,probes,nearGenes,percentage=0.2,permu.size=10000,
                     permu.dir=NULL, Pe=0.001,dir.out="./",diffExp=FALSE,cores=NULL,
                     portion = 0.3, label=NULL,save=TRUE){
  ## check data
  if(!all(probes %in% rownames(mee@meth))) 
    stop("Probes option should be subset of rownames of methylation matrix.")
  if(is.character(nearGenes)){
    newenv <- new.env()
    load(nearGenes, envir=newenv)
    nearGenes <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }else if(!is.list(nearGenes)){
    stop("nearGene option must be a list containing output of GetNearGenes function 
         or path of rda file containing output of GetNearGenes function.")
  }
  #get raw pvalue
  ##I need to modify that if there is all NA. stop the process.
  if(requireNamespace("parallel", quietly=TRUE) && requireNamespace("snow", quietly=TRUE)) {
	  if(!is.null(cores)){
		  if(cores > parallel::detectCores()) cores <- parallel::detectCores()/2
		  cl <- snow::makeCluster(cores,type = "SOCK")
		  Probe.gene<-parallel::parSapplyLB(cl,probes,Stat.nonpara,Meths= getMeth(mee,probe=probes), 
											NearGenes=nearGenes,K=portion,Top=percentage,
											Exps=getExp(mee),simplify = FALSE)
		  parallel::stopCluster(cl)
	  }else{
	    Probe.gene<-sapply(probes,Stat.nonpara,Meths=getMeth(mee,probe=probes),
	                       NearGenes=nearGenes,K=portion,Top=percentage,Exps=mee@exp,
	                       simplify = FALSE)
	  }
  } else {
    Probe.gene<-sapply(probes,Stat.nonpara,Meths=getMeth(mee,probe=probes),
                       NearGenes=nearGenes,K=portion,Top=percentage,Exps=mee@exp,
                       simplify = FALSE)
  }
  
  Probe.gene <- do.call(rbind,Probe.gene)

  Probe.gene <- Probe.gene[!is.na(Probe.gene$Raw.p),]
  #   Probe.gene$logRaw.p <- -log10(Probe.gene$Raw.p)
  GeneID <- unique(Probe.gene[!is.na(Probe.gene$Raw.p),"GeneID"])
  # get permutation
  permu <- get.permu(mee,geneID=GeneID, percentage=percentage, rm.probes=probes, 
                     permu.size=permu.size, portion = portion,
                     permu.dir=permu.dir,cores=cores)
  #get empirical p-value
  message("Calculate empirical P value.\n")
  Probe.gene.Pe <- Get.Pvalue.p(Probe.gene,permu)
  Probe.gene.Pe <- Probe.gene.Pe[order(Probe.gene.Pe$Raw.p),]
  if(save) write.csv(Probe.gene.Pe, file=sprintf("%s/getPair.%s.all.pairs.statistic.csv",
                                        dir.out, label),row.names=FALSE)
  selected <- Probe.gene.Pe[Probe.gene.Pe$Pe < Pe & !is.na(Probe.gene.Pe$Pe),]
  
  if(diffExp){
    ## calculate differential expression between two groups.
    Exp <- getExp(mee, geneID = unique(selected$GeneID))
    TN <- getSample(mee,cols = "TN")
    out <- lapply(split(Exp,rownames(Exp)),
                    function(x, TN){test <- t.test(x~TN)
                                    U.test <- wilcox.test(x~TN)
                                    out <- data.frame("log2FC_TvsN" = test$estimate[2]-test$estimate[1],
                                                      "TN.diff.pvalue"=test$p.value)
                                    return(out)}, TN=TN)
    out <- do.call(rbind, out)
    out$GeneID <- rownames(out)
    add <- out[match(selected$GeneID, out$GeneID),c("log2FC_TvsN","TN.diff.pvalue")]
    selected <- cbind(selected, add)                                                         
  }
  if(save) write.csv(selected, file=sprintf("%s/getPair.%s.pairs.significant.csv",dir.out, label),
            row.names=FALSE)
  invisible(gc())
  return(selected)
}

### permutation
#permu.size can be all which mean all the usable probes.
#' get.permu
#' @param mee A MEE.data object containing at least meth, exp, probeInfo, geneInfo.
#' @param geneID A vector lists the genes' ID.
#' @param rm.probes A vector lists the probes name.
#' @param cores A interger which defines number of core to be used in parallel process.
#'  Default is NULL: don't use parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of 
#' samples used to link probes to genes. Default is 0.2.
#' @param permu.size A number specify the times of permuation. Default is 1000.
#' @param permu.dir A path where the output of permuation will be. 
#' @param portion A number specify the cut point for methylated and unmethylated.
#' Default is 0.3.
#' @return Permutations
#' @export 
#' @examples
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' permu <-get.permu(mee=mee,geneID=rownames(getExp(mee)),
#'                   rm.probes=c("cg00329272","cg10097755"),
#'                   permu.size=5)
get.permu <- function(mee, geneID, percentage=0.2, rm.probes=NULL ,portion=0.3,
                      permu.size=10000, permu.dir=NULL,cores=NULL){
  set.seed(200)
  ## get usable probes
  binary.m <- rowMeans(Binary(mee@meth,portion),na.rm = TRUE)
  usable.probes <- names(binary.m[binary.m <0.95 & binary.m > 0.05 & !is.na(binary.m)])
  usable.probes <- usable.probes[!usable.probes %in% rm.probes]
  if(length(usable.probes) < permu.size) 
    stop(sprintf("There is no enough usable probes to perform %s time permutation, 
                 set a smaller permu.size.",permu.size))
  if(!is.numeric(permu.size)) permu.size <- length(usable.probes) 
  probes.permu <- sample(usable.probes, size = permu.size, replace = FALSE)
  
  if(is.null(permu.dir)){
	  permu.meth <- getMeth(mee,probe=probes.permu)
	  if(requireNamespace("parallel", quietly=TRUE) && requireNamespace("snow", quietly=TRUE)) {
		  if(!is.null(cores)){
			  if(cores > parallel::detectCores()) cores <- parallel::detectCores()/2
			  suppressWarnings(cl <- snow::makeCluster(cores,type = "SOCK"))
			  permu<-parallel::parSapplyLB(cl,probes.permu,Stat.nonpara.permu,Meths=permu.meth,
										   Gene=unique(as.character(getGeneInfo(mee)$GENEID)),
										   Top=percentage,Exps=getExp(mee), 
										   simplify = FALSE)
			  parallel::stopCluster(cl)
		  } else {
			  permu<-sapply(probes.permu,Stat.nonpara.permu,Meths=permu.meth,
							Gene=unique(as.character(getGeneInfo(mee)$GENEID)),
							Top=percentage,Exps=getExp(mee),
							simplify=FALSE)
		  }
	  } else {
		  permu<-sapply(probes.permu,Stat.nonpara.permu,Meths=permu.meth,
						Gene=unique(as.character(getGeneInfo(mee)$GENEID)),
						Top=percentage,Exps=getExp(mee),
						simplify=FALSE)
	  }
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
		  permu.meth <- getMeth(mee,probe=tmp.probes)
		  if(requireNamespace("parallel", quietly=TRUE) && requireNamespace("snow", quietly=TRUE)) {
			  if(!is.null(cores)){
				  if(cores > parallel::detectCores()) cores <- parallel::detectCores()/2
				  suppressWarnings(cl <- snow::makeCluster(cores,type = "SOCK"))
				  permu<-parallel::parSapplyLB(cl,tmp.probes,Stat.nonpara.permu,Meths=permu.meth,
											   Gene=unique(as.character(getGeneInfo(mee)$GENEID)),
											   Top=percentage,Exps=getExp(mee), permu.dir=permu.dir,
											   simplify = FALSE)
				  parallel::stopCluster(cl)
			  } else {
				  permu<-sapply(tmp.probes,Stat.nonpara.permu,Meths=permu.meth,
								Gene=unique(as.character(getGeneInfo(mee)$GENEID)),
								Top=percentage,Exps=getExp(mee),permu.dir=permu.dir,
								simplify=FALSE)
			  }
		  } else {
			  permu<-sapply(tmp.probes,Stat.nonpara.permu,Meths=permu.meth,
							Gene=unique(as.character(getGeneInfo(mee)$GENEID)),
							Top=percentage,Exps=getExp(mee),permu.dir=permu.dir,
							simplify=FALSE)
		  }
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
#' @import S4Vectors IRanges GenomicRanges
#' @return A data frame contains genes whose expression significantly anti-correlated 
#' with promoter methylation.
#' @export
promoterMeth <- function(mee,sig.pvalue=0.01,percentage=0.2,save=TRUE){
  TSS_2K <- promoters(getGeneInfo(mee), upstream = 100, downstream = 700)
  probes <- getProbeInfo(mee)
  overlap <- findOverlaps(probes, TSS_2K)
  df <- data.frame(Probe=as.character(probes$name[queryHits(overlap)]), 
                   GeneID=TSS_2K$GENEID[subjectHits(overlap)], stringsAsFactors=FALSE)
  if(nrow(df)==0){
    out <- data.frame(GeneID=c(), Symbol=c(), Raw.p= c())
  }else{
    df <- unique(df)
    ProbeInTSS <- split(df$Probe,df$GeneID)
    
    ## calculate average methylation of promoter
    Gene.promoter <- lapply(ProbeInTSS, 
                            function(x, METH){meth <- METH[x,]
                                              if(length(x)>1){
                                                meth <- colMeans(meth,na.rm=TRUE)
                                              }  
                                              return(meth)},   
                            METH=getMeth(mee))
    Gene.promoter <- do.call(rbind, Gene.promoter)
    ## make fake NearGene 
    Fake <- data.frame(Symbol = getSymbol(mee, geneID = rownames(Gene.promoter)),
                       GeneID = rownames(Gene.promoter),
                       Distance= 1,
                       Side = 1, stringsAsFactors=FALSE)
    Fake <- split(Fake, Fake$GeneID)
    out <- lapply(rownames(Gene.promoter),Stat.nonpara, NearGenes=Fake,K=0.3,Top=0.2,
                  Meths=Gene.promoter,Exps=getExp(mee))
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
#' enriched.motif <- get.enriched.motif(probes=probes,background.probes = bg,
#' min.incidence=2, label="hypo")
get.enriched.motif <- function(probes.motif, probes, background.probes,
							   lower.OR=1.1,min.incidence=10, dir.out="./",
							   label=NULL,save=TRUE){
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
			newenv <- new.env()
			load(sprintf("%s/probeInfo_feature_distal.rda",dir.out), envir=newenv)
      background.probes <- get(ls(newenv)[1],envir=newenv) 
      background.probes <- as.character(background.probes$name)
      # The data is in the one and only variable
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
  Summary <- data.frame(motif = colnames(probes.TF), NumOfProbes= probes.TF.num,
                        OR=sub.enrich.TF, lowerOR=sub.enrich.TF.lower, 
                        upperOR=sub.enrich.TF.upper)
  Summary <- Summary[order(Summary$lowerOR, decreasing = TRUE),]
  if(save) write.csv(Summary, file= sprintf("%s/getMotif.%s.motif.enrichment.csv",
                                   dir.out,label))
  
  ## enriched motif and probes
  en.motifs <- names(sub.enrich.TF.lower[sub.enrich.TF.lower > lower.OR &
                                           !sub.enrich.TF.lower %in% "Inf" & 
                                           probes.TF.num > min.incidence])
  message(sprintf("%s motifs are enriched.",length(en.motifs)))
  enriched.motif <- sapply(en.motifs, 
                           function(x, probes.TF)
                             {names(probes.TF[probes.TF[,x]==1,x])},
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

#' get.TFs
#' @param mee A MEE.data object containing at least meth, exp, probeInfo, geneInfo. 
#' @param enriched.motif Can be either a list containing output of 
#' get.enriched.motif function or path of rda file containing output of get.enriched.motif function.
#' @param TFs Can be either a data.frame containing TF GeneID and Symbol 
#' or path of csv file containing TF GeneID and Symbol. If missing, 
#' human TF list will be used. For detail information, refer reference paper.
#' @param motif.relavent.TFs Can be either a list containing motif (list name) 
#' and relavent TF (content of list) or path of rda file containing a list 
#' containing motif (list name) and relavent TF (content of list). If missing,
#'  human TF list will be used. For detail information, refer reference paper.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of 
#' samples used to link probes to genes. Default is 0.2.
#' @param cores A interger which defines number of core to be used in parallel 
#' process. Default is NULL: don't use parallel process.
#' @param dir.out A path specify the directory for outputs. Default is current directory
#' @param label A character labels the outputs.
#' @return Potential responsible TFs will be reported.
#' @export 
#' @examples
#' load(system.file("extdata","mee.example.rda",package = "ELMER"))
#' enriched.motif <- list("TP53"= c("cg00329272", "cg10097755", "cg08928189",
#'                                  "cg17153775", "cg21156590", "cg19749688", "cg12590404",
#'                                  "cg24517858", "cg00329272", "cg09010107", "cg15386853",
#'                                  "cg10097755", "cg09247779", "cg09181054"))
#'TF <- get.TFs(mee, enriched.motif, 
#'               TFs=data.frame(GeneID=c("ID7157","ID8626","ID7161"),
#'               Symbol=c("TP53","TP63","TP73"), 
#'               stringsAsFactors = FALSE),
#'               label="hypo")
get.TFs <- function(mee, enriched.motif, TFs, motif.relavent.TFs,
                    percentage=0.2,dir.out="./",label=NULL,cores=NULL,save=TRUE){
  if(missing(enriched.motif)){
    stop("enriched.motif is empty.")
  }else if(is.character(enriched.motif)){
    newenv <- new.env()
    load(enriched.motif, envir=newenv)
    enriched.motif <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }else if(!is.list(enriched.motif)){
    stop("enriched.motif option should be a list object.")
  }
  
  if(missing(TFs)){
    newenv <- new.env()
    data("human.TF",package = "ELMER.data",envir=newenv)
    TFs <- get(ls(newenv)[1],envir=newenv) 
    if(all(grepl("ID", rownames(getExp(mee)) )))
      TFs$GeneID <- paste0("ID",TFs$GeneID)
    TFs <- TFs[TFs$GeneID %in% rownames(getExp(mee)),]
  }else if(is.character(TFs)){
    TFs <- read.csv(TFs, stringsAsFactors=FALSE)
  }
  
  if(missing(motif.relavent.TFs)){
    newenv <- new.env()
    data("motif.relavent.TFs",package = "ELMER.data",envir=newenv)
    motif.relavent.TFs <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }else if(is.character(motif.relavent.TFs)){
    newenv <- new.env()
    load(motif.relavent.TFs, envir=newenv)
    motif.relavent.TFs <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }
  
  motif.meth <- lapply(enriched.motif, 
                       function(x,meth){if(length(x)<2)
                         { return(meth[x,])
                       }else{
                         return(colMeans(meth[x,],na.rm = TRUE))
                       }}, meth = getMeth(mee,probe=unique(unlist(enriched.motif))) )
                                                                   
  motif.meth <- do.call(rbind, motif.meth)

  if(requireNamespace("parallel", quietly=TRUE) && requireNamespace("snow", quietly=TRUE)) {
	  if(!is.null(cores)){
		  if(cores > parallel::detectCores()) cores <- parallel::detectCores()/2
		  cl <- snow::makeCluster(cores,type = "SOCK")
		  TF.meth.cor<-parallel::parSapplyLB(cl,names(enriched.motif),
											 Stat.nonpara.permu,Meths=motif.meth,Gene=TFs$GeneID,
											 Top=percentage,Exps=getExp(mee), simplify=FALSE)
		  parallel::stopCluster(cl)
	  }else{
	    TF.meth.cor<-sapply(names(enriched.motif),Stat.nonpara.permu,Meths=motif.meth,
	                        Gene=TFs$GeneID,Top=percentage,Exps=getExp(mee), simplify=FALSE) 
	  }
  } else {
    TF.meth.cor<-sapply(names(enriched.motif),Stat.nonpara.permu,Meths=motif.meth,
                        Gene=TFs$GeneID,Top=percentage,Exps=getExp(mee), simplify=FALSE) 
  }
  TF.meth.cor <- lapply(TF.meth.cor, function(x){return(x$Raw.p)})
  TF.meth.cor <- do.call(cbind,TF.meth.cor)
  ## check row and col names
  rownames(TF.meth.cor) <- TFs$Symbol
  cor.summary <- sapply(colnames(TF.meth.cor), 
                        function(x, TF.meth.cor, motif.relavent.TFs)
                          { cor <- sort(TF.meth.cor[,x])
                            top <- names(cor[1:floor(0.05*nrow(TF.meth.cor))])
                            potential.TF <- top[top %in% motif.relavent.TFs[[x]]]
                            out <- data.frame("motif"=x,"top potential TF"= potential.TF[1],
                                              "potential TFs"= paste(potential.TF, collapse = ";"),
                                              "top_5percent"= paste(top,collapse = ";"))},                                         
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
