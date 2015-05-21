
\docType{methods}
\name{getGeneInfo}
\alias{getGeneInfo}
\alias{getGeneInfo,ANY-method}
\title{getGeneInfo to extract geneInfo slot from MEE.data or Pair object.}
\description{
getGeneInfo is a function to easily extract geneInfo out of a MEE.data or Pair object.
By specifying geneID or symbol, geneInfo for the defined genes (geneID or symbol) will be
extracted out of MEE.data or Pair object. When range is specified, the geneInfo falling into
the range will be extracted.
}

\usage{
getGeneInfo(object, geneID, symbol, range)

\S4method{getGeneInfo}{ANY}(object, geneID, symbol, range)
}
\arguments{
\item{object}{MEE.data or Pair object}

\item{geneID}{A vector of genes' id. When specified, only these gene
coordinates will be output.}

\item{symbol}{A vector of genes' symbols . When specified, only these gene
coordinates will be output.}

\item{range}{A GRanges object. When specified, only the geneInfo locating within
these loci will be output.}
}
\value{
Gene annotation information such as gene id, symbol and coordinates.
}

\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}

\examples{
geneInfo <- txs()
mee <- fetch.mee(geneInfo=geneInfo)
Genes <- getGeneInfo(mee, geneID = "55811")
Genes <- getGeneInfo(mee, symbol ="ADCY10")
Genes <- getGeneInfo(mee, range = GRanges(seqnames="chr1", ranges=IRanges(1000000,1600000)))
}

