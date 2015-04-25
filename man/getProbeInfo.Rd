
\docType{methods}
\name{getProbeInfo}
\alias{getProbeInfo}
\alias{getProbeInfo,ANY-method}
\title{getProbeInfo to extract probeInfo slot from MEE.data or Pair object.}
\description{
getProbeInfo is a function to easily extract probeInfo out of a MEE.data or Pair object.
By specifying probe, probeInfo for the defined set of probes will be extracted out of 
MEE.data or Pair object. Option chr will restrain the output probeInfo to certain chromosomes.
When range is specified, the probeInfo falling into the range will be extracted.
}

\usage{
getProbeInfo(object, chr, probe, range)

\S4method{getProbeInfo}{ANY}(object, chr, probe, range)
}
\arguments{
\item{object}{MEE.data or Pair object}

\item{chr}{A vector of chromosome such chr1, chr2. When specified, only the
probeInfo locating on these chromosome will be output.}

\item{probe}{A vector of probes' name. When specified, only these probes'
coordinate will be output.}

\item{range}{A GRanges object. When specified, only probeInfo locating
within these loci will be output.}
}
\value{
Probe information such as names, coordinates.
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}
\examples{
probeInfo <- GRanges(seqnames = c("chr1","chr1","chr3"),
ranges = IRanges(start = c(1,6,20),end = c(2,7,21)),
name=c("cg1","cg2","cg3"))
mee <- fetch.mee(probeInfo=probeInfo)
Probes <- getProbeInfo(mee,chr="chr1") # get probes which locate on the chr1
Probes <- getProbeInfo(mee, probe = "cg1") # get certain probes information
Probes <- getProbeInfo(mee, range = GRanges(seqnames="chr1", ranges=IRanges(5,20)))
}

