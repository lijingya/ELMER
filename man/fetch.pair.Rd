\name{fetch.pair}
\alias{fetch.pair}
\title{
 fetch.pair to generate Pair class object.
}
\description{
  fetch.pair is a funtion to take in enhancer-gene linkage prediction information,
  probe information and gene annotation generating a Pair class object, which is the 
  input for plotting functions. Options (pair, probeInfo, geneInfo) can
  take in R object or read files by specifying file paths. 
}
\usage{
fetch.pair(pair, probeInfo, geneInfo)
}

\arguments{
  \item{pair}{
 A data.frame (R object) or a path of XX.csv file containing pair information such as
 output of function \code{\link{get.pair}}.
}
  \item{probeInfo}{
 A GRnage object or a path of XX.rda file which only contains a GRange of probe information.
}
  \item{geneInfo}{
 A GRnage object or path of XX.rda file which only contains a GRange of gene 
 information such as Coordinates, GENEID and SYMBOL. 
}
}


\value{
A Pair class object containing 3 slots. Detail see \code{\link{Pair-class}}
}

\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}



\seealso{
\code{\link{schematic.plot}}
}
\examples{
df <- data.frame(Probe=c("cg19403323","cg12213388","cg26607897"),
GeneID =c("ID255928","ID84451","ID55811"),
Symbol =c("SYT14","KIAA1804","ADCY10"),
Pe=c(0.003322259,0.003322259,0.003322259))
geneInfo <- txs()
## input can be a path
pair <- fetch.pair(pair = df, geneInfo=geneInfo)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
