# Stat for computing hypo- and hypermethylated probes
##ID: the rownames of Exps
##TN: the two group factors
##Exps: the data matrix. Rows is the ID and cols are sample which divided into two groups.
##cutoff: filter out the small difference.
#' Stat for computing hypo- and hypermethylated probes
#' @param Probe A character which is name of probes on the array.
#' @param TN A vector of characters either 'Tumor' or 'Normal' which represent samples from Tumor or Normal.
#' @param test A test method used for the statistic.
#' @param Tumor.per A number determines the percentage of Tumor samples will be used in the test.
#' @param Normal.per A number determines the percentage of Normal samples will be used in the test.
#' @param hyper A character either NULL,FALSE, TRUE. NULL indicate two way test to identify hypo/hypermethylated probes. FALSE indicate to identify hypomethylated probes. TRUE indicates to identify hypermethylated probes.
#' @return test results.
.Stat <- function(Probe,TN,Meth,test=t.test,Tumor.per=NULL,Normal.per=NULL,hyper=NULL){
  out <- c()
  Meth <- Meth[Probe,]
  if(!is.na(Tumor.per)){
    if(Top.m){
      tumor.tmp <- sort(Meth[TN %in% "Tumor"],decreasing=T)
      normal.tmp <- sort(Meth[TN %in% "Normal"],decreasing=T)
    }else{
      tumor.tmp <- sort(Meth[TN %in% "Tumor"])
      normal.tmp <- sort(Meth[TN %in% "Normal"])
    }
    if(round(length(normal.tmp)*Normal.per)< 5){
      if(length(normal.tmp) < 5) {
        Normal.number <- length(normal.tmp)
      }else{
        Normal.number <- 5
        message(sprintf("%s percentage of normal sample is less than 5. Set number of normal samples as 5",Normal.per))
      }
    }else{
      Normal.number <- round(length(normal.tmp)*Normal.per)
    }
    tumor.tmp <- tumor.tmp[1:round(length(tumor.tmp)*Tumor.per)]
    normal.tmp <- normal.tmp[1:Normal.number]
    Meth <- c(normal.tmp,tumor.tmp)
    TN <- c(rep("Normal",length(normal.tmp)),rep("Tumor",length(tumor.tmp)))
  }
  ##this is to remove the situation that the normal or tumor are all NA (only one is value)
  Meth_split <- split(Meth,TN)
  Meth_split <- unlist(lapply(Meth_split,function(x){!is.na(sd(x,na.rm=T))}))
  
  if(sd(Meth,na.rm=T)>0 & all(Meth_split)){
    if(!is.null(hyper)){
      alternative <- ifelse(hyper,"less","greater")
    }else{
      alternative <- "two.sided"
    }
    df <- data.frame(Meth=Meth,TN=factor(TN))
    TT <- test(Meth~TN,df,alternative=alternative)
    MeanDiff <- TT$estimate[2]-TT$estimate[1]
    PP <- TT$p.value
    out <- rbind(out,c(Probe,PP,MeanDiff))
  }else{
    out <- rbind(out,c(Probe,NA,NA))
  }
  return(out)
}





#---different statistic test----------------------------
##Probes must be one probes
##Gene can be one gene or multiple gene.
#' U test (non parameter test) for permutation. This is one probe vs multiple gene which is good for computing permutation for each probe.
#' @param Probe A character of name of Probe in array.
#' @param Gene A vector of gene ID.
#' @param K A number determines the methylated groups and unmethylated groups.
#' @param Top A number determines the percentage of top methylated/unmethylated samples.
#' @param Meths A matrix contains methylation for each probe (row) and each sample (column).
#' @param Exps A matrix contains Expression for each gene (row) and each sample (column).
#' @return U test results
   
.Stat.nonpara.permu <- function(Probe,Gene,K,Top=NULL,Meths=Meths,Exps=Exps){
  if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
  Exp <- as.matrix(Exps[Gene,])
  Meth <- Meths[Probe,]
  Meth_B <- Binary(Meth,Break=K)
  test.p <- c()
  if(sum(Exp)==0){
    test.p <- NA
  }else{
    unmethyPercent <- sum(Meth_B==0,na.rm=T)/length(Meth_B)
    methyPercent <- sum(Meth_B==1,na.rm=T)/length(Meth_B)
    if(unmethyPercent < 0.05 | methyPercent < 0.05){
      test.p <- NA
    }else{
      unmethy <- order(Meth)[1:round(length(Meth)*Top)] 
      methy <- order(Meth,decreasing=T)[1:round(length(Meth)*Top)] 
      Fa <- factor(rep(NA,length(Meth)),levels=c(-1,1))
      Fa[unmethy] <- -1
      Fa[methy] <- 1
      if(length(Gene)<2){
        test.p <- wilcox.test(Exp~Fa,alternative = "greater",exact=F)$p.value
      }else{
        cl <- makeCluster(13,"SOCK")
        #        test.p <- unlist(mclapply(1:nrow(Exp),function(x,Factor) { exp = Exp[x,]
        #                                                           wilcox.test(exp~Factor,alternative = "greater",exact=F)$p.value},Factor=Fa,mc.cores=10))
        test.p <- unlist(parSapplyLB(cl,1:nrow(Exp),function(x,Factor) { exp = Exp[x,]
                                                                         wilcox.test(exp~Factor,alternative = "greater",exact=F)$p.value},Factor=Fa,simplify=F))
        stopCluster(cl)
        
      }  
    }
    if(length(Gene) < 2){
      out <- c(Probe,Gene,test.p)
    }else{
      out <- cbind(Probe=rep(Probe,length(Gene)),Gene,test.p)
    }
    
  } 
  return(out)
}

#' U test (non parameter test) for permutation. This is one probe vs nearby gene which is good for computing each probes for nearby genes.
#' @param Probe A character of name of Probe in array.
#' @param NearGenes A list of nearby gene for each probe which is output of GetNearGenes function.
#' @param K A number determines the methylated groups and unmethylated groups.
#' @param Top A number determines the percentage of top methylated/unmethylated samples.
#' @param Meths A matrix contains methylation for each probe (row) and each sample (column).
#' @param Exps A matrix contains Expression for each gene (row) and each sample (column).
#' @return U test results
.Stat.nonpara <- function(Probe,NearGenes,K,Top=NULL,Meths=Meths,Exps=Exps){
  source("/export/uec-gs1/laird/users/lijingya/software/scripts/R/Heatmap.Func.R")
  if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
  Gene <- NearGenes[[Probe]][,2]
  Exp <- as.matrix(Exps[Gene,])
  Meth <- Meths[Probe,]
  Meth_B <- Binary(Meth,Break=K)
  test.p <- c()
  if(sum(Exp)==0){
    test.p <- NA
  }else{
    unmethyPercent <- sum(Meth_B==0,na.rm=T)/length(Meth_B)
    methyPercent <- sum(Meth_B==1,na.rm=T)/length(Meth_B)
    if(unmethyPercent < 0.05 | methyPercent < 0.05){
      test.p <- NA
    }else{
      unmethy <- order(Meth)[1:round(length(Meth)*Top)] 
      methy <- order(Meth,decreasing=T)[1:round(length(Meth)*Top)] 
      Fa <- factor(rep(NA,length(Meth)),levels=c(-1,1))
      Fa[unmethy] <- -1
      Fa[methy] <- 1
      test.p <- apply(Exp,1,function(x,Factor) {wilcox.test(x~Factor,alternative = "greater",exact=F)$p.value},Factor=Fa)
    }
    out <- cbind(Probe=rep(Probe,length(Gene)),Gene,test.p)
  } 
  return(out)
}

#' Calculate empirical Pvalue
#' @param zscore.Matrix A data.frame of raw pvalue from U test. Output from .Stat.nonpara
#' @param permu data frame of permutation. Output from .Stat.nonpara.permu
#' @return A data frame with empirical Pvalue.
Get.Pvalue.p <- function(zscore.Matrix,permu){
  .Pvalue <- function(x,permu,target.Matrix){
    zscore <- target.Matrix[x,4]
    Gene <- target.Matrix[x,"Gene"]
    if(is.na(zscore)){
      out <- NA
      #print("NA")
    }else{
      out <- (sum(permu[Gene,]  > zscore | permu[Gene,] == zscore,na.rm=T)+1)/(sum(!is.na(permu[Gene,])) + 1)
    }
    #  else if(zscore <= 0){
    #    out <- sum(permu[Gene,] < zscore | permu[Gene,] == zscore, na.rm=T)/length(permu[Gene,])   
    return(out)
  }
  cl <- makeCluster(5,"SOCK")
  Pvalue <- parSapplyLB(cl,1:nrow(zscore.Matrix),.Pvalue,target.Matrix=zscore.Matrix,permu=permu,simplify = F)
  stopCluster(cl)
  Pvalue <- unlist(Pvalue)
  Output <- cbind(zscore.Matrix,Pvalue)
  return(Output)
}
