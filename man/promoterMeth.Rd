
\name{promoterMeth}
\alias{promoterMeth}
\title{promoterMeth to calculate associations of gene expression with DNA methylation
at promoter regions.}
\usage{
promoterMeth(mee, sig.pvalue = 0.01, percentage = 0.2, save = TRUE)
}
\arguments{
\item{mee}{A MEE.data object must contains four components: meth, exp, probeInfo, geneInfo.}

\item{sig.pvalue}{A number specifies significant cutoff for gene silenced by promoter
methylation. Default is 0.01. P value is raw P value without adjustment.}

\item{percentage}{A number ranges from 0 to 1 specifying the percentage of
samples of control and experimental groups used to link promoter DNA methylation to genes.
Default is 0.2.}

\item{save}{A logic. If it is true, the result will be saved.}
}

\value{
A data frame contains genes whose expression significantly anti-correlated
with promoter methylation.
}
\description{
promoterMeth is a function to calculate associations of gene expression with DNA methylation
at promoter regions.
}

\references{
Yao L, Shen H, Laird PW, Farnham PJ,Berman BP: Inferring Regulatory Element 
Landscapes and Transcription Factor Networks from Cancer Methylomes. in revision of
Genome Biology
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}

\examples{
load(system.file("extdata","mee.example.rda",package = "ELMER"))
Gene.promoter <- promoterMeth(mee) 
}