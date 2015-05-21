## get differential methylated probes-------------------------
#' Stat.diff.meth
#' @param probe A charactor specify probe name
#' @param meths A matrix contain DNA methylation data.
#' @param TN A vector of category of samples.
#' @param test A function specify which statistic test will be used.
#' @param percentage A number specify the percentage of normal and tumor 
#' samples used in the test.
#' @param Top.m A logic. If to identify hypomethylated probe Top.m should be FALSE. 
#' hypermethylated probe is TRUE.
#' @return Statistic test results to identify differentially methylated probes.
Stat.diff.meth <- function(probe,meths,TN,test=t.test,percentage=0.2,Top.m=NULL){
  meth <- meths[probe,]
  if(Top.m){
    tumor.tmp <- sort(meth[TN %in% "Experiment"],decreasing=TRUE)
    normal.tmp <- sort(meth[TN %in% "Control"],decreasing=TRUE)
  }else{
    tumor.tmp <- sort(meth[TN %in% "Experiment"])
    normal.tmp <- sort(meth[TN %in% "Control"])
  }
  if(round(length(normal.tmp)*percentage)< 5){
    if(length(normal.tmp) < 5) {
      Normal.number <- length(normal.tmp)
    }else{
      Normal.number <- 5
    }
  }else{
    Normal.number <- round(length(normal.tmp)*percentage)
  }
  tumor.tmp <- tumor.tmp[1:round(length(tumor.tmp)*percentage)]
  normal.tmp <- normal.tmp[1:Normal.number]
  meth <- c(normal.tmp,tumor.tmp)
  TN <- c(rep("Control",length(normal.tmp)),rep("Experiment",length(tumor.tmp)))
  
  ##this is to remove the situation that the normal or tumor are all NA (only one is value)
  meth_split <- split(meth,TN)
  meth_split <- unlist(lapply(meth_split,function(x){!is.na(sd(x,na.rm=TRUE))}))
  
  if(sd(meth,na.rm=TRUE)>0 & all(meth_split)){
    if(!is.na(Top.m)){
      alternative <- ifelse(Top.m,"less","greater")
    }else{
      alternative <- "two.sided"
    }
    df <- data.frame(meth=meth,TN=factor(TN))
    TT <- test(meth~TN,df,alternative=alternative)
    MeanDiff <- TT$estimate[2]-TT$estimate[1]
    PP <- TT$p.value
    out <- data.frame(probe=probe,PP=PP,MeanDiff=MeanDiff, stringsAsFactors = FALSE)
  }else{
    out <- data.frame(probe=probe,PP=NA,MeanDiff=NA,stringsAsFactors = FALSE)
  }
  return(out)
}

#'Stat.nonpara.permu
#' @param Probe A character of name of Probe in array.
#' @param Gene A vector of gene ID.
#' @param Top A number determines the percentage of top methylated/unmethylated samples.
#' @param Meths A matrix contains methylation for each probe (row) and each sample (column).
#' @param Exps A matrix contains Expression for each gene (row) and each sample (column).
#' @param permu.dir A path to store permuation data.
#' @return U test results
Stat.nonpara.permu <- function(Probe,Gene,Top=0.2,Meths=Meths,Exps=Exps,permu.dir=NULL){
  if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
  Exp <- Exps[Gene,]
  if(is.vector(Meths)){
    Meth <- Meths
  }else{
    Meth <- Meths[Probe,]
  }
  unmethy <- order(Meth)[1:round(length(Meth)*Top)] 
  methy <- order(Meth,decreasing=TRUE)[1:round(length(Meth)*Top)] 
  Fa <- factor(rep(NA,length(Meth)),levels=c(-1,1))
  Fa[unmethy] <- -1
  Fa[methy] <- 1
  Exp <- Exp[,!is.na(Fa)]
  Fa <- Fa[!is.na(Fa)]
  test.p <- unlist(lapply(splitmatrix(Exp),
                          function(x,Factor) 
                            {wilcox.test(x[Factor %in% -1],x[Factor %in% 1],alternative = "greater",exact=FALSE)$p.value},
							Factor=Fa))
  out <- data.frame(GeneID=Gene,
                    Raw.p=test.p[match(Gene, names(test.p))], 
                    stringsAsFactors = FALSE) 
  if(is.null(permu.dir)){
    return(out)
  }else{
    write.table(out,file=sprintf("%s/%s",permu.dir,Probe),
                quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  }
}


#' U test (non parameter test) for permutation. This is one probe vs nearby gene 
#' which is good for computing each probes for nearby genes.
#' @param Probe A character of name of Probe in array.
#' @param NearGenes A list of nearby gene for each probe which is output of GetNearGenes function.
#' @param K A number determines the methylated groups and unmethylated groups.
#' @param Top A number determines the percentage of top methylated/unmethylated samples.
#' @param Meths A matrix contains methylation for each probe (row) and each sample (column).
#' @param Exps A matrix contains Expression for each gene (row) and each sample (column).
#' @return U test results
Stat.nonpara <- function(Probe,NearGenes,K,Top=NULL,Meths=Meths,Exps=Exps){
  if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
  Gene <- NearGenes[[Probe]][,2]
  Exp <- Exps[Gene,]
  Meth <- Meths[Probe,]
  Meth_B <- mean(Binary(Meth,Break=K),na.rm = TRUE)
  if( Meth_B >0.95 | Meth_B < 0.05 ){
    test.p <- NA
  }else{
    unmethy <- order(Meth)[1:round(length(Meth)*Top)] 
    methy <- order(Meth,decreasing=TRUE)[1:round(length(Meth)*Top)] 
    Fa <- factor(rep(NA,length(Meth)),levels=c(-1,1))
    Fa[unmethy] <- -1
    Fa[methy] <- 1
    if(!is.vector(Exp)){
      Exp <- Exp[,!is.na(Fa)]
      Fa <- Fa[!is.na(Fa)]
      test.p <- unlist(lapply(splitmatrix(Exp),
                              function(x,Factor) 
                              {wilcox.test(x[Factor %in% -1],x[Factor %in% 1],
                                           alternative = "greater",
                                           exact=FALSE)$p.value},
                              Factor=Fa))
    }else{
      Exp <- Exp[!is.na(Fa)]
      Fa <- Fa[!is.na(Fa)]
      test.p <- wilcox.test(Exp[Fa %in% -1],Exp[Fa %in% 1],
                                           alternative = "greater",
                                           exact=FALSE)$p.value
    }
  }
  
  if(length(Gene)==1){
    out <- data.frame(Probe=rep(Probe,length(Gene)),
                      GeneID=Gene,Symbol=NearGenes[[Probe]]$Symbol, 
                      Distance=NearGenes[[Probe]]$Distance, 
                      Sides=NearGenes[[Probe]]$Side,
                      Raw.p=test.p, 
                      stringsAsFactors = FALSE)
  }else{
    out <- data.frame(Probe=rep(Probe,length(Gene)),
                      GeneID=Gene,Symbol=NearGenes[[Probe]]$Symbol, 
                      Distance=NearGenes[[Probe]]$Distance, 
                      Sides=NearGenes[[Probe]]$Side,
                      Raw.p=test.p[match(Gene, names(test.p))], 
                      stringsAsFactors = FALSE)
  }
  
  return(out)
}


#' Calculate empirical Pvalue
#' @param U.matrix A data.frame of raw pvalue from U test. Output from .Stat.nonpara
#' @param permu data frame of permutation. Output from .Stat.nonpara.permu
#' @return A data frame with empirical Pvalue.
Get.Pvalue.p <- function(U.matrix,permu){
  .Pvalue <- function(x,permu){
    Raw.p <- as.numeric(x["Raw.p"])
    Gene <- as.character(x["GeneID"])
    if(is.na(Raw.p)){
      out <- NA
    }else{
      out <- (sum(permu[as.character(Gene),]  < Raw.p | 
                    permu[Gene,] == Raw.p,na.rm=TRUE)+1)/(sum(!is.na(permu[Gene,])) + 1)
    } 
    return(out)
  }
  Pvalue <- unlist(apply(U.matrix,1,.Pvalue,permu=permu))
  U.matrix$Pe <- Pvalue
  return(U.matrix)
}
