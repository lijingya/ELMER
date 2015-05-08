
\name{getGeneID}
\alias{getGeneID}
\title{getGeneID to report gene id from symbol}
\description{
getGeneID is a function to report the gene ids from gene symbols.
}
\usage{
getGeneID(mee, symbol)
}
\arguments{
\item{mee}{A MEE.data or Pair object.}

\item{symbol}{A vector of characters which are gene symbols }
}
\value{
The gene ID for these gene symbols
}

\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}

\examples{
geneInfo <- txs()
## input can be a path
pair <- fetch.pair(geneInfo=geneInfo)
getGeneID(pair, symbol="KIAA1804")
}

