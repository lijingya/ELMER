\name{get.feature.probe}
\alias{get.feature.probe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.feature.probe to select probes within promoter regions or distal regions.
}
\description{
get.feature.probe is a function to select the probes falling into distal feature regions or
promoter regions.
}
\usage{
get.feature.probe(feature, TSS, TSS.range = list(upstream = 2000, downstream = 2000), promoter = FALSE, rm.chr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{feature}{
A GRange object containing biofeature coordinates such as enhancer coordinates. 
Default is comprehensive genomic enhancer regions from REMC and FANTOM5 which is
Union.enhancer data in \pkg{ELMER.data}. feature option is only usable when 
promoter option is FALSE.
}
  \item{TSS}{
A GRange object contains the transcription start sites. When promoter is FALSE, Union.TSS
in \pkg{ELMER.data} will be used for default. When promoter is TRUE, UCSC gene TSS will
be used as default (see detail). User can specify their own preference TSS annotation. 
}
  \item{TSS.range}{
A list specify how to define promoter regions. Default is upstream =2000bp 
and downstream=2000bp.
}
  \item{promoter}{
A logical. If it is TRUE, function will ouput the promoter probes.
If FALSE, function will ouput the distal probes overlaping with features. The
default is FALSE.
}
  \item{rm.chr}{
A vector of chromosome. Once specified, the probes on these chromosome will be removed
such as chrX chrY or chrM
}
}
\details{
  TSS:
    In order to get real distal probes, we use more comprehensive annotated TSS by both 
    GENCODE and UCSC. However, to get probes within promoter regions need more
    accurate annotated TSS such as UCSC. Therefore, there are different settings for
    promoter and distal probe selection. But user can specify their own favorable
    TSS annotation. Then there won't be any difference between promoter and distal
    probe selection.
  
}
\value{
A GRanges object contains the coordinate of probes which locate within promoter regions
or distal feature regions such as union enhancer from REMC and FANTOM5.
}

\references{
Yao L, Shen H, Laird PW, Farnham PJ,Berman BP: Inferring Regulatory Element 
Landscapes and Transcription Factor Networks from Cancer Methylomes. in revision of
Genome Biology
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}


\examples{
  # get distal enhancer probe
\dontrun{
Probe <- get.feature.probe()
}
# get promoter probes
\dontrun{
Probe <- get.feature.probe(promoter=TRUE)
}
# get distal enhancer probe remove chrX chrY
Probe2 <- get.feature.probe(rm.chr=c("chrX", "chrY"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
