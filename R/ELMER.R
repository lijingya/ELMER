#' ELMER is designed to use DNA methylation and gene expression from a
#' large number of samples to infere regulatory element landscape and transcription
#' factor network in primary tissue.
#' @docType package
#' @name ELMER
NULL

#' A MultiAssayExperiment containing
#' DNA methylation data: 101 probes from platform 450K 
#' Gene Expression data: 1026 genes
#' for 234 samples from TCGA-LUSC.
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name data
#' @format A MultiAssayExperiment for 234 Samples (8 normal samples, 226 Primary solid tumor)
NULL

#' A matrix containing DNA methylation beta-values from TCGA
#' DNA methylation data: 1728 probes
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name Meth
#' @format A MultiAssayExperiment for 268 Samples  and 1728 probes
NULL

#' A matrix containing gene expression data from TCGA
#' Gene Expression data: 3842 genes
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name GeneExp
#' @format A gene expression matrix for 234 Samples and 3842 genes
NULL


#' A matrix containing the relashionship between the DNA methylation level at the probes
#' for a given motif and the gene expression of TF.
#' This data is used in the examples of ELMER package
#' @docType data
#' @keywords internal
#' @name TF.meth.cor
#' @format A matrix with 1968 rows (TFs) and 12 columns (motifs)
NULL

