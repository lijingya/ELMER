\docType{methods}
\name{getExp}
\alias{getExp}
\alias{getExp,MEE.data-method}
\title{getExp to extract RNA expression slot from a MEE.data class object.}
\description{
getExp is a function to easily extract RNA expression data out of a MEE.data object.
By specifying geneID or ID for samples, a matrix of RNA expression data for defined
genes and samples will be extract out of MEE.data object.
}
\usage{
getExp(object, geneID, ID)

\S4method{getExp}{MEE.data}(object, geneID, ID)
}
\arguments{
\item{object}{MEE.data object}

\item{geneID}{A vector of genes' id. When specified, only these gene expressions 
will be output.}

\item{ID}{A vector of sample ID. When specified, gene expression only for
these samples will be output.}
}
\value{
A matrix of gene expression values. 
}

\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}

\examples{
exp <- matrix(data=c(101:110),ncol=5,dimnames=list(c("gene1","gene2"),paste0("sample",1:5)))
mee <- fetch.mee(exp=exp)
Exp <- getExp(mee, geneID = "gene1") ## get gene expression for certain genes
Exp <- getExp(mee, ID = c("sample1","sample5")) ## get gene expression for certain samples
}

