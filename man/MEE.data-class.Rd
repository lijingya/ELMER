
\docType{class}
\name{MEE.data-class}
\alias{MEE.data-class}
\title{MEE.data
An S4 class contains 5 slots: methylation, expression, sample information, probe
information and gene information. MEE.data class are the main input for main functions.}
\description{
MEE.data
An S4 class contains 5 slots: methylation, expression, sample information, probe
information and gene information. MEE.data class are the main input for main functions.
}
\section{Slots}{

\describe{
\item{\code{meth}}{A matrix of DNA methylation. Each row is one probe and each
column is one sample}

\item{\code{exp}}{A matrix of expression. Each row is one gene and each
column is one sample}

\item{\code{sample}}{A data.frame contains sample information}

\item{\code{probeInfo}}{A GRange object contains probe information}

\item{\code{geneInfo}}{A GRange object contains gene information}
}}

\seealso{
  \code{\link{get.diff.meth}}, \code{\link{get.pair}}, \code{\link{get.TFs}}, \code{\link{get.permu}}, 
  \code{\link{getMeth}},\code{\link{getExp}}, \code{\link{getGeneInfo}}, \code{\link{getProbeInfo}},
  \code{\link{getSample}},\code{\link{fetch.mee}},\code{\link{scatter.plot}}
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}
