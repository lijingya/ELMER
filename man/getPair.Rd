
\docType{methods}
\name{getPair}
\alias{getPair}
\alias{getPair,Pair-method}
\title{getPair to extract pairInfo slot from Pair object.}
\description{
getPair is a function to easily extract pairInfo out of a Pair object.
By specifying geneID or probe, geneInfo for the defined genes (geneID ) and probes (probe)
will be extracted out of Pair object. 
}
\usage{
getPair(object, probe, geneID)

\S4method{getPair}{Pair}(object, probe, geneID)
}
\arguments{
\item{object}{Pair object}

\item{probe}{A vector of probes' name. When specified, only the pair containing
these probes will be output.}

\item{geneID}{A vector of genes' id. When specified, only the pair containing
these genes will be output.}
}
\value{
Pair information such as empirical P values, probe and gene ID.
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}
\examples{
df <- data.frame(Probe=c("cg19403323","cg12213388","cg26607897"),
GeneID =c("ID255928","ID84451","ID55811"),
Symbol =c("SYT14","KIAA1804","ADCY10"),
Pe=c(0.003322259,0.003322259,0.003322259))
pair <- fetch.pair(pair = df)
Pairs <- getPair(pair, probe = "cg19403323") # get pair information for a probe
Pairs <- getPair(pair, geneID = "ID55811") # get pair information for a gene
}

