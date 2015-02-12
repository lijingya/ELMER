##get.distal.en
#' get.feature.probe
#' @param probe A GRange object containing probes coordinate information. Default is Illumina-methyl-450K probes coordinates.
#' @param distal A logical. If FALSE, function will output the all probes overlaping with features. If TRUE, function will ouput the distal probes overlaping with features.
#' @param feature A GRange object containing biofeature coordinate such as enhancer coordinates. Default is comprehensive genomic enhancer regions from REMC and FANTOM5.
#' @param TSS A GRange object containing the transcription start site. Default is UCSC gene TSS.
#' @param TSS.range A list specify how to define promoter regions. Default is upstream =2000bp and downstream=2000bp.
#' @param rm.chr A vector of chromosome need to be remove from probes such as chrX chrY or chrM
#' @return A GRange object containing probes that satisfy selecting critiria.
#' @export 

get.feature.probe <- function(probe,distal=TRUE,feature,TSS,TSS.range=list(upstream=2000,downstream=2000),rm.chr=NULL){
  if(missing(probe)){
    warning("Default probes coordinates are for HM450K DNA methylation array")
    probe <- ReadBed(system.file("extdata","Illumina-methyl-450K-manifest.hg19.bed",package = "ELMER"))
  }
  if(!is.null(rm.chr)) probe <- probe[!as.character(seqnames(probe)) %in% rm.chr]
  if(distal){
    if(missing(TSS)){
      newenv <- new.env()
      load(system.file("extdata","UCSC_gene_hg19.rda",package = "ELMER"),envir=newenv)
      txs <- get(ls(newenv)[1],envir=newenv)
      TSS <- promoters(txs,upstream = TSS.range[["upstream"]], downstream = TSS.range[["downstream"]])    
    }else{
      TSS <- promoters(TSS,upstream = TSS.range[["upstream"]], downstream = TSS.range[["downstream"]])
    }
    probe <- probe[setdiff(1:length(probe),unique(queryHits(findOverlaps(probe,TSS))))]
  }
  
  
  if(missing(feature)){
    feature <- ReadBed(system.file("extdata","Union_strong_enhancer_REMC_FANTOM.bed",package = "ELMER"))
  }else if(!is(feature,"GRange")){              ## feature can be specified as file path or GRange
    feature <- ReadBed(feature)
  }
  probe <- probe[unique(queryHits(findOverlaps(probe,feature)))]  
  return(probe)
}

## get differential methylated probes-------------------------
## TCGA pipe don't specify dir.out
#' get.diff.meth
#' @param mee A MEE.data object containing at least meth and probeInfo.
#' @param diff.dir A character showing differential methylation dirction. It can be "hypo" which is only selecting hypomethylated probes; "hyper" which is only selecting hypermethylated probes; "both" which select both hypomethyalted and hypermethylated probes.
#' @param cores A interger which defines number of core to be used in parallel process. Default is NULL: don't use parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples used to identify the differential methylation. Default is 0.2.
#' @param pvalue A number specify the significant Pvalue cutoff for significant hypo/hyper-methylated probes. Default is current directory.
#' @param sig.dif A number specify the significant methylation difference cutoff for significant hypo/hyper-methylated probes. Default is 0.3.
#' @param dir.out A path specify the directory for outputs. Default is is current directory.
#' @return Statistics for all probes and significant hypo or hyper-methylated probes.
#' @export 
get.diff.meth <- function(mee,diff.dir="both",cores=NULL,percentage=0.2, pvalue=0.01, sig.dif=0.3, dir.out="./"){
  if(nrow(mee@meth)==0) stop("Cannot identify differential DNA methylation region without DNA methylation data.")
  if(nrow(getSample(mee))==0){
    stop("Sample information data to do differential analysis.")
  }else if(is.null(getSample(mee,cols="TN"))){
    stop("\"TN\" should be specified, labeling two group of sample for comparison.")
  }else if(length(table(getSample(mee,cols="TN")))<2){
    stop("\"TN\" should have at 2 distinct group labels for comparison.")
  }
  if("both" == diff.dir) diff.dir <- c("hypo","hyper")
  
  result <- list()
  if(!is.null(cores)){
    if(require(snow)){
      require(parallel)
      if(cores > detectCores()) cores <- detectCores()/2
      cl <- makeCluster(cores,type = "SOCK")
    }
  }
  if("hyper" %in% diff.dir){
    if(!is.null(cores)){
      out <- parSapplyLB(cl,rownames(mee@meth),t.Stat,percentage=percentage,meth=mee@meth,TN=getSample(mee,cols="TN"),Top.m=TRUE,simplify =F)
    }else{
      out <- sapply(rownames(mee@meth),t.Stat,percentage=percentage,meth=mee@meth,TN=getSample(mee,cols="TN"),Top.m=TRUE,simplify =F)
    }
    out <- do.call(rbind,out)
    out <- as.data.frame(out,stringsAsFactors = F)
    out$adjust.p <- p.adjust(as.numeric(out[,2]),method="BH")
    colnames(out) <- c("probe","pvalue","tumorMinNormal","adjust.p")
    write.csv(out,file=sprintf("%s/hyper.probes.csv",dir.out), row.names=FALSE)
    write.csv(out[out$adjust.p < pvalue & abs(out$tumorMinNormal)>sig.dif,],file=sprintf("%s/hyper.probes_significant.csv",dir.out), row.names=FALSE)
    result[["hyper"]] <- out[out$adjust.p < pvalue & abs(out$tumorMinNormal)>sig.dif,]
  }
  if("hypo" %in% diff.dir){
    if(!is.null(cores)){
      out <- parSapplyLB(cl,rownames(mee@meth),t.Stat,percentage=percentage,meth=mee@meth,TN=getSample(mee,cols="TN"),Top.m=FALSE,simplify =F)
    }else{
      out <- sapply(rownames(mee@meth),t.Stat,percentage=percentage,meth=mee@meth,TN=getSample(mee,cols="TN"),Top.m=FALSE,simplify =F)
    }
    out <- do.call(rbind,out)
    out <- as.data.frame(out,stringsAsFactors = F)
    out$adjust.p <- p.adjust(as.numeric(out[,2]),method="BH")
    colnames(out) <- c("probe","pvalue","tumorMinNormal","adjust.p")
    write.csv(out,file=sprintf("%s/hypo.probes.csv",dir.out), row.names=FALSE)
    write.csv(out[out$adjust.p < pvalue & abs(out$tumorMinNormal)>sig.dif,],file=sprintf("%s/hypo.probes_significant.csv",dir.out), row.names=FALSE)
    result[["hypo"]] <- out[out$adjust.p < pvalue & abs(out$tumorMinNormal)>sig.dif,]
  }
  
  if(!is.null(cores)){
    stopCluster(cl)
  }
  
   
  return(result)  
}


### pairing function--------------------
#nearGenes : can be either a list containing output of GetNearGenes function or path of rda file containing output of GetNearGenes function.
## TCGA pipe don't specify dir.out
#' get pair
#' @param mee A MEE.data object containing at least meth, exp, probeInfo, geneInfo.
#' @param probes A vector lists the probes' name that need to be linked to genes.
#' @param nearGenes Can be either a list containing output of GetNearGenes function or path of rda file containing output of GetNearGenes function.
#' @param cores A interger which defines number of core to be used in parallel process. Default is NULL: don't use parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples used to link probes to genes. Default is 0.2.
#' @param permu.size A number specify the times of permuation. Default is 1000.
#' @param permu.dir A path where the output of permuation will be. 
#' @param Pe A number specify the empircal pvalue cutoff for defining signficant pairs.
#' @param dir.out A path specify the directory for outputs. Default is current directory
#' @param label A character labels the outputs.
#' @return Statistics for all pairs and significant pairs
#' @export 
get.pair <- function(mee,probes,nearGenes,percentage=0.2,permu.size=1000,permu.dir=NULL, Pe=0.01,dir.out="./",cores=NULL,label=NULL){
  ## check data
  browser()
  if(!all(probes %in% rownames(mee@meth))) stop("Probes option should be subset of rownames of methylation matrix.")
  if(is.character(nearGenes)){
    newenv <- new.env()
    load(nearGenes, envir=newenv)
    nearGenes <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }else if(!is.list(nearGenes)){
    stop("nearGene option must be a list containing output of GetNearGenes function or path of rda file containing output of GetNearGenes function.")
  }
  if(is.null(permu.dir)) permu.dir <- paste0(dir.out,"/","permu")
  #get raw pvalue
  ##I need to modify that if there is all NA. stop the process.
  if(!is.null(cores)){
    if(require(snow)){
      require(parallel)
      if(cores > detectCores()) cores <- detectCores()/2
      cl <- makeCluster(cores,type = "SOCK")
    }
    Probe.gene<-parSapplyLB(cl,probes,Stat.nonpara,NearGenes=nearGenes,K=0.3,Top=percentage,Meths=getMeth(mee,probe=probes),
                            Exps=getExp(mee,GeneID=unique(unlist(lapply(nearGenes,function(x){paste0("ID",x[,"GeneID"])})))),simplify = F)
    stopCluster(cl)
  }else{
    Probe.gene<-sapply(probes,Stat.nonpara,NearGenes=nearGenes,K=0.3,Top=percentage,Meths=mee@meth,Exps=mee@exp,simplify = F)
  }
  
  Probe.gene <- do.call(rbind,Probe.gene)
  Probe.gene <- Probe.gene[!is.na(Probe.gene$Raw.p),]
#   Probe.gene$logRaw.p <- -log10(Probe.gene$Raw.p)
  GeneID <- unique(Probe.gene[!is.na(Probe.gene$Raw.p),"GeneID"])
  # get permutation
  permu <- get.permu(mee,GeneID=GeneID, percentage=percentage, rm.probes=probes, permu.size=100, permu.dir=permu.dir,cores=6)
  #get empirical p-value
  message("Calculate empirical P value.\n")
  Probe.gene.Pe <- Get.Pvalue.p(Probe.gene,permu)
  Probe.gene.Pe <- Probe.gene.Pe[order(Probe.gene.Pe$Raw.p),]
  write.csv(Probe.gene.Pe, file=sprintf("%s/%s.all.pairs.statistic.csv",dir.out, label),row.names=F)
  selected <- Probe.gene.Pe[Probe.gene.Pe$Pe < Pe & !is.na(Probe.gene.Pe$Pe),]
  write.csv(selected, file=sprintf("%s/%s.significant.pairs.csv",dir.out, label),row.names=F)
}

### permutation
#permu.size can be all which mean all the usable probes.
#' get.permu
#' @param mee A MEE.data object containing at least meth, exp, probeInfo, geneInfo.
#' @param GeneID A vector lists the genes' ID.
#' @param rm.probes A vector lists the probes name.
#' @param cores A interger which defines number of core to be used in parallel process. Default is NULL: don't use parallel process.
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples used to link probes to genes. Default is 0.2.
#' @param permu.size A number specify the times of permuation. Default is 1000.
#' @param permu.dir A path where the output of permuation will be. 
#' @return permutation
#' @export 
get.permu <- function(mee, GeneID=GeneID, percentage=0.2, rm.probes=NULL ,permu.size=1000, permu.dir=NULL,cores=NULL){
  set.seed(200)
  ## get usable probes
  binary.m <- rowMeans(Binary(mee@meth,0.3),na.rm = T)
  usable.probes <- names(binary.m[binary.m <0.95 & binary.m > 0.05])
  usable.probes <- usable.probes[!usable.probes %in% rm.probes]
  if(length(usable.probes) < permu.size) stop(sprintf("There is no enough usable probes to perform %s time permutation, set a smaller permu.size.",permu.size))
  usable.probes <- sample(usable.probes,size = permu.size,replace = F)
  if(!is.numeric(permu.size)) permu.size <- length(usable.probes) 
  probes.permu <- sample(usable.probes, size = permu.size, replace = F)
  
  if(!is.null(cores)){
    if(require(snow)){
      require(parallel)
      if(cores > detectCores()) cores <- detectCores()/2
      suppressWarnings(cl <- makeCluster(cores,type = "SOCK"))
    }
  }
  ## if file already there don't need to calculate.
  if(!file.exists(permu.dir)){
    dir.create(permu.dir,recursive = T)
  }
  if(!all(probes.permu %in% dir(permu.dir))){
    tmp.probes <- probes.permu[!probes.permu %in% dir(permu.dir)]
    if(!is.null(cores)){
      permu<-parSapplyLB(cl,tmp.probes,Stat.nonpara.permu,Gene=as.character(getGeneInfo(mee)$GENEID),Top=percentage,Meths=getMeth(mee,probe=tmp.probes),Exps=getExp(mee), permu.dir=permu.dir,simplify = F)
    }else{
      permu<-sapply(tmp.probes,Stat.nonpara.permu,Gene=as.character(getGeneInfo(mee)$GENEID),Top=percentage,Meths=getMeth(mee,probe=tmp.probes),Exps=getExp(mee),permu.dir=permu.dir,simplify = F)
    }
  }
  permu.p <- paste0(permu.dir,"/",probes.permu)
  if(!is.null(cores)){
    permu <- parSapplyLB(cl,permu.p,function(x,GeneID){ tmp <- read.table(x,stringsAsFactors=F)
                                                        tmp <- tmp[match(GeneID,tmp[,1]),2]},GeneID=GeneID,simplify=F)
  }else{
    permu <- sapply(permu.p,function(x,GeneID){ tmp <- read.table(x,stringsAsFactors=F)
                                                tmp <- tmp[match(GeneID,tmp[,1]),2]},GeneID=GeneID,simplify=F)
  }
  stopCluster(cl)
  permu <- do.call(cbind,permu)
  rownames(permu) <- GeneID
  colnames(permu) <- probes.permu
  return(permu)
}



#' get.enriched.motif
#' @param probes.motif A matrix contains motifs occurrence within probes regions. 
#' @param probes A vector lists the probes' names in which motif enrichment will be calculated.
#' @param background.probes A vector list of probes' names which are considered as background for motif.enrichment calculation.
#' @param lower.OR A number specify the lower boundary of Odds ratio which defines the significant enriched motif.
#' @param min.incidence A non-negative integer specify the minimum incidence of motif in the given probes set.
#' @param dir.out A path specify the directory for outputs. Default is current directory
#' @param label A character labels the outputs.
#' @return A list contains enriched motifs with the probes regions harboring the motif.
#' @export 
get.enriched.motif <- function(probes.motif, probes, background.probes,lower.OR=1.1,min.incidence=10, dir.out="./",label=NULL){
  browser()
  if(missing(probes.motif)){
    probes.motif <- system.file("extdata","Probesall.TF.matrix.200bp.rda",package = "ELMER")
    newenv <- new.env()
    load(probes.motif, envir=newenv)
    all.probes.TF <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  } 
  ## here need to be add motif search part.
  if(missing(probes)) stop("probes option should be specified.")
  if(missing(background.probes)){
    if(file.exists(sprintf("%s/probeInfo_feature.rda",dir.out))){
      newenv <- new.env()
      load(sprintf("%s/probeInfo_feature.rda",dir.out), envir=newenv)
      background.probes <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
      background.probes <- as.character(background.probes$name)
    }else{
      background.probes <- colnames(probes.motif)
    }
  }
  Probes.TF.percent <- colMeans(Probes.TF)
  ## load probes for enriched motif ----------------------------------------------
  probes.TF <- all.probes.TF[probes,]
  probes.TF.num <- colSums(probes.TF, na.rm=T)
  sub.enrich.TF <- colMeans(probes.TF)*(1-Probes.TF.percent)/Probes.TF.percent/(1-colMeans(probes.TF))
  SE <- sqrt(1/colSums(probes.TF) + 1/(nrow(probes.TF)-colSums(probes.TF)) + 1/colSums(probes.TF)+ 1/(nrow(probes.TF)-colSums(probes.TF)))
  sub.enrich.TF.lower <- exp(log(sub.enrich.TF)-1.96*SE)
  sub.enrich.TF.upper <- exp(log(sub.enrich.TF)+1.96*SE)
  ## summary
  Summary <- data.frame(motif = colnames(probes.TF), NumOfProbes= probes.TF.num, OR=sub.enrich.TF, lowerOR=sub.enrich.TF.lower, upperOR=sub.enrich.TF.upper)
  Summary <- Summary[order(Summary$lowerOR, decreasing = T),]
  write.csv(Summary, file= sprintf("%s/%s.motif.enrichment.csv"))
  
  ## enriched motif and probes
  en.motifs <- names(sub.enrich.TF.lower[sub.enrich.TF.lower > lower.OR & !sub.enrich.TF.lower %in% "Inf" & probes.TF.num > min.incidence])
  cat(sprintf("%s motifs are enriched.",length(en.motifs)))
  enriched.motif <- sapply(en.motifs, function(x, probes.TF){names(probes.TF[probes.TF[,x]==1,x])}, probes.TF=probes.TF)
  save(enriched.motif, file= sprintf("%s/%s.enriched.motifs.rda"))
  return(enriched.motif)
  
  ## add information to siginificant pairs
  if(file.exists(sprintf("%s/%s.significant.pairs.csv",dir.out, label))){
    sig.Pairs <- read.csv(sprintf("%s/%s.significant.pairs.csv",dir.out, label), stringsAsFactors=F)
    if(all(sig.Pairs$Probes %in% rownames(Probes.TF))){
      motif.Info <- sapply(sig.Pairs$Probes, function(x, Probes.TF,en.motifs){TFs <- names(Probes.TF[x,Probes.TF[x,]==1])
                                                                              non.en.motif <- paste(setdiff(TFs,en.motifs),collapse = ";")
                                                                              en.motif <- paste(intersect(TFs,en.motifs), collapse = ";")
                                                                              out <- data.frame(non_enriched_motifs=non.en.motif, enriched_motifs=en.motif, stringsAsFactors = F)
                                                                              return(out)},Probes.TF=Probes.TF, en.motifs=en.motifs)
      motif.Info <- do.call(rbind,motif.Info)
      sig.Pairs <- cbind(sig.Pairs, motif.Info)
      write.csv(selected, file=sprintf("%s/%s.significant.pairs.withmotif.csv",dir.out, label),row.names=F)
    }
  }
}


get.TFs <- function(mee, enriched.motif, TFs, motif.relavent.TFs){
  if(missing(enriched.motif)) stop("enriched.motif is empty.")
  if(missing(TFs)) TFs <- read.csv(system.file("extdata","human.TF.list.csv",package = "ELMER"), stringsAsFactors=F)
  if(missing(motif.relavent.TFs)){
    newenv <- new.env()
    load(system.file("extdata","motif.relavent.TFs.human.rda",package = "ELMER"), envir=newenv)
    motif.relavent.TFs <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }
  if(is.character(enriched.motif)){
    newenv <- new.env()
    load(enriched.motif, envir=newenv)
    enriched.motif <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }
  if(is.character(TFs) && length(TFs)==1){
    TFs <- read.csv(TFs, stringsAsFactors=F)
  }
  if(is.character(motif.relavent.TFs)){
    newenv <- new.env()
    load(motif.relavent.TFs, envir=newenv)
    motif.relavent.TFs <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }
  
}