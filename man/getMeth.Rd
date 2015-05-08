
\docType{methods}
\name{getMeth}
\alias{getMeth}
\alias{getMeth,MEE.data-method}
\title{getMeth to extract DNA methylation slot from a MEE.data class object.}
\description{
getMeth is a function to easily extract DNA methylation data out of a MEE.data object.
By specifying probe or ID for samples, a matrix of DNA methylation values for defined
probes and samples will be extracted out of MEE.data object.
}
\usage{
getMeth(object, probe, ID)

\S4method{getMeth}{MEE.data}(object, probe, ID)
}
\arguments{
\item{object}{MEE.data object}

\item{probe}{A vector of probes' name. When specified, DNA methylation only
for these probes will be output.}

\item{ID}{A vector of sample ID. When specified, DNA methylation only for
these samples will be output.}
}
\value{
A matrix of DNA methylation values.
}

\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}

\examples{
meth <- matrix(data=c(1:20),ncol=5,dimnames=list(paste0("probe",1:4),paste0("sample",1:5)))
mee <- fetch.mee(meth=meth)
Meth <- getMeth(mee,probe = "probe1")
Meth <- getMeth(mee, ID = c("sample1","sample2"))
}

