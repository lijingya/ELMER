\name{GetNearGenes}
\alias{GetNearGenes}
\title{GetNearGenes to collect nearby genes for one locus.}
\description{
GetNearGenes is a function to collect equal number of gene on each side of one locus.
}
\usage{
GetNearGenes(geneAnnot = NULL, TRange = NULL, geneNum = 20,
  cores = NULL)
}
\arguments{
\item{geneAnnot}{A GRange object contains coordinates of promoter for
human genome.}

\item{TRange}{A GRange object contains coordinates of a list of target loci.}

\item{geneNum}{A number determines how many gene will be collected totally. Then
the number devided by 2 is the number of genes collected from each side of targets 
(number shoule be even) Default to 20.}

\item{cores}{A interger which defines the number of cores to be used in parallel process.
Default is NULL: no parallel process.}
}
\value{
A data frame of nearby genes' information: genes' IDs, genes' symbols,
distance with target and side to which the gene locate to the target.
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
geneAnnot <- txs(TSS=list(upstream=0, downstream=0))
probe <- GRanges(seqnames = c("chr1","chr2"),
range=IRanges(start = c(16058489,236417627), end= c(16058489,236417627)),
name= c("cg18108049","cg17125141"))
NearbyGenes <- GetNearGenes(geneNum=20,geneAnnot=geneAnnot,TRange=probe)
}

