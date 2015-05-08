
\name{txs}
\alias{txs}
\title{txs to fetch USCS gene annotation (transcripts level) from 
Bioconductor package Homo.sapians. If upstream and downstream are specified in TSS list, promoter
regions of USCS gene will be generated.}
\usage{
txs(TSS = list(upstream = NULL, downstream = NULL))
}
\arguments{
\item{TSS}{A list contains upstream and downstream like TSS=list(upstream, downstream).
When upstream and downstream is specified, coordinates of promoter regions with
gene annotation will be generated.}
}
\value{
UCSC gene annotation if TSS is not specified. Coordinates
of UCSC gene promoter regions if TSS is specified.
}
\description{
txs is a function to fetch USCS gene annotation (transcripts level) from 
Bioconductor package Homo.sapians. If upstream and downstream are specified in 
TSS list, promoter regions of USCS gene will be generated.
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}
\examples{
# get UCSC gene annotation (transcripts level)
txs <- txs()
# get coordinate of all UCSC promoter regions +/-1000bp of TSSs
\dontrun{
txs <- txs(TSS=list(upstream=1000, downstream=1000))
}
}

