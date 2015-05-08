\name{fetch.mee}
\alias{fetch.mee}

\title{
  fetch.mee to generate MEE.data class object.
}
\description{
  fetch.mee is a funtion to take in DNA methylation, RNA expression, sample information
  probe information and gene annotation generating a MEE.data class object, which is the 
  input for main functions. Options (meth, exp, sample, probeInfo, geneInfo) can
  take in R object or read files by specifying file paths. When TCGA is specified, sample
  information will be automatically generated such as Control/Experiment labels. 
}
\usage{
fetch.mee(meth, exp, sample, probeInfo, geneInfo, probes = NULL, genes = NULL, TCGA = FALSE)
}

\arguments{
  \item{meth}{
 A matrix (R object) or path of XX.rda file which only stores a matrix of
 DNA methylation data.
}
  \item{exp}{
 A matrix (R object) or path of XX.rda file which only stores a matrix of
 expression data.
}
  \item{sample}{
 A data frame (R object) or path of XX.rda file which only contains sample
 information in data frame format.
}
  \item{probeInfo}{
 A GRnage object or path of XX.rda file which only contains a GRange of probe 
 information.
}
  \item{geneInfo}{
 A GRnage object or path of XX.rda file which only contains a GRange of gene 
 information such as Coordinates, GENEID and SYMBOL. 
}
  \item{probes}{
 A vector lists name of probes. If probes are specified, the DNA methylation matrix and 
 probeInfo in MEE.data object will be restrained to this list of probes.
}
  \item{genes}{
 A vector lists gene ids. If gene are specified, the methylation and probeInfo 
 in output MEE.data object will be restrained this list of genes.
}
  \item{TCGA}{
 A logical. FALSE indicates that data is not from TCGA (FALSE is default).
 TRUE indicates data is from TCGA and sample section will automatically filled in.
}
}

\note{
 Options (meth, exp, sample, probeInfo, geneInfo) don't need to be all specified.
 User can input one or more according the needs of the function. Such as get.meth.diff
 only need meth, sample and probeInfo.
}
\details{
 Use path to load in data will help to reduce memory usage.
}
\value{
A MEE.data object containing 5 slots. Detail see \code{\link{MEE.data-class}}
}

\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}


\examples{
meth <- matrix(data=c(1:20),ncol=5,dimnames=list(paste0("probe",1:4),paste0("sample",1:5)))
exp <- matrix(data=c(101:110),ncol=5,dimnames=list(c("gene1","gene2"),paste0("sample",1:5)))
mee <- fetch.mee(meth=meth, exp=exp)
## only fetch probe 1 and 3
mee <- fetch.mee(meth=meth, exp=exp, probes=c("probe1","probe3"))
## only fetch gene 1
mee <- fetch.mee(meth=meth, exp=exp, genes="gene1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
