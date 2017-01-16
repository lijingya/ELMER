## get differential methylated probes-------------------------
#' Stat.diff.meth
#' @param probe A charactor specify probe name
#' @param meth A matrix contain DNA methylation data.
#' @param TN A vector of category of samples.
#' @param test A function specify which statistic test will be used.
#' @param percentage A number specify the percentage of normal and tumor 
#' samples used in the test.
#' @param Top.m A logic. If to identify hypomethylated probe Top.m should be FALSE. 
#' hypermethylated probe is TRUE.
#' @return Statistic test results to identify differentially methylated probes.
Stat.diff.meth <- function(probe,
                           meth,
                           groups,
                           group1,
                           group2,
                           test=t.test,
                           percentage=0.2,
                           Top.m=NULL){
  group1.tmp <- sort(meth[groups %in% group1],decreasing = Top.m)
  group2.tmp <- sort(meth[groups %in% group2],decreasing = Top.m)
  group1.nb <- ifelse(round(length(group1.tmp) * percentage) < 5, min(5,length(group1.tmp)), round(length(group1.tmp) * percentage))
  group2.nb <- ifelse(round(length(group2.tmp) * percentage) < 5, min(5,length(group2.tmp)), round(length(group2.tmp) * percentage))
  
  group1.tmp <- group1.tmp[1:group1.nb]
  group2.tmp <- group2.tmp[1:group2.nb]
  
  if(sd(meth,na.rm=TRUE)>0 & !all(is.na(group1.tmp)) & !all(is.na(group2.tmp))){
    if(!is.na(Top.m)){
      alternative <- ifelse(Top.m,"greater","less")
    } else {
      alternative <- "two.sided"
    }
    # If hyper (top. TRUE alternative greater) group 1 > group 2
    # If hypo  (top. FALSE alternative greater) group 1 < group 2
    TT <- test(x = group1.tmp, y = group2.tmp, alternative = alternative, conf.int = TRUE)
    
    if(length(TT$estimate) == 2) {
      MeanDiff <- TT$estimate[1]-TT$estimate[2]
    } else {
      MeanDiff <- TT$estimate
    }
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
#' @return U test results
Stat.nonpara.permu <- function(Probe,
                               Gene,
                               Top=0.2,
                               Meths=Meths,
                               Exps=Exps){
  idx <- order(Meths)
  nb <- round(length(Meths)*Top)
  unmethy <- head(idx, n = nb) 
  methy <- tail(idx, n = nb) 
  
  test.p <- unlist(lapply(splitmatrix(Exps),
                          function(x) {
                            wilcox.test(x[unmethy],x[methy],alternative = "greater",exact=FALSE)$p.value
                          }))
  
  test.p <- data.frame(GeneID=Gene,
                       Raw.p=test.p[match(Gene, names(test.p))], 
                       stringsAsFactors = FALSE) 
  return(test.p)
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
Stat.nonpara <- function(Probe,
                         NearGenes,
                         K,
                         Top=NULL,
                         Meths=Meths,
                         Exps=Exps){
  if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
  Gene <- NearGenes[[Probe]][,2]
  Exp <- Exps[Gene,]
  Meth <- Meths[Probe,]
  Meth_B <- mean(Meth > K, na.rm = TRUE)
  if( Meth_B >0.95 | Meth_B < 0.05 ){
    test.p <- NA
  }else{
    idx <- order(Meth)
    nb <- round(length(Meth)*Top)
    unmethy <- head(idx, n = nb) 
    methy <- tail(idx, n = nb) 
    # Here we will test if the Expression of the unmethylated group is higher than the exptression of the methylated group
    if(!is.vector(Exp)){
      Exps <- Exps[,c(unmethy,methy)]
      test.p <- unlist(lapply(splitmatrix(Exp),
                              function(x,Factor) 
                              {wilcox.test(x[unmethy],x[methy],
                                           alternative = "greater",
                                           exact=FALSE)$p.value},
                              Factor=Fa))
    }else{
      test.p <- wilcox.test(Exps[unmethy],Exps[methy],
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
    } else {
      #       num( Pp ≤ Pr) + 1
      # Pe = ---------------------
      #            x + 1
      # Pp = pvalue probe (Raw.p)
      # Pr = pvalue random probe (permu matrix)
      out <- (sum(permu[as.character(Gene),]  <= Raw.p, na.rm=TRUE) + 1) / (sum(!is.na(permu[Gene,])) + 1)
    } 
    return(out)
  }
  message("Calculate empirical P value.\n")
  Pvalue <- unlist(apply(U.matrix,1,.Pvalue,permu=permu))
  U.matrix$Pe <- Pvalue
  return(U.matrix)
}
