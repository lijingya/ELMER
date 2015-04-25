
\name{getSymbol}
\alias{getSymbol}
\title{getSymbol to report gene symbol from id}
\usage{
getSymbol(mee, geneID)
}
\arguments{
\item{mee}{A MEE.data or Pair object.}

\item{geneID}{A character which is the geneID}
}
\value{
The gene symbol for input genes.
}
\description{
getSymbol is a function to report the gene symbols from gene IDs.
}

\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}

\examples{
geneInfo <- txs()
## input can be a path
pair <- fetch.pair(geneInfo=geneInfo)
getSymbol(pair, geneID="84451")
}

