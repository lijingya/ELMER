
\name{scatter.plot}
\alias{scatter.plot}
\title{scatter.plot to plot scatter plots between gene expression and DNA methylation.}
\usage{
scatter.plot(mee, byPair = list(probe = c(), gene = c()),
  byProbe = list(probe = c(), geneNum = 20), byTF = list(TF = c(), probe =
  c()), category = NULL, dir.out = "./", save = TRUE, ...)
}
\arguments{
\item{mee}{A MEE.data object includes DNA methylation data, expression data,
probeInfo and geneInfo.}

\item{byPair}{A list: byPair =list(probe=c(),gene=c()); probe contains a vector
of probes' name and gene contains a vector of gene IDs. The length of probe
should be the same with length of gene.}

\item{byProbe}{A list byProbe =list(probe=c(), geneNum=20); probe contains
a vector of probes' name and geneNum specify the number of gene near the probes
will ploted. 20 is default for geneNum.}

\item{byTF}{A list byTF =list(TF=c(), probe=c()); TF contains a vector of TF's
symbol and probe contains the a vector of probes' name. }

\item{category}{A vector labels subtype of samples or a character which is the
column name in the sampleInfo in the MEE.data object. Once specified, samples
will label different color. The color can be customized by using color.value.}

\item{dir.out}{A path specifies the directory to which the figures will be saved.
Current directory is default.}

\item{save}{A logic. If true, figure will be saved to dir.out.}

\item{...}{color.value, lm_line in scatter function}
}
\value{
Scatter plots.
}
\description{
scatter.plot is a function to plot various scatter plots between gene expression and 
DNA methylation. When byPair is specified, scatter plot for individual probe-gene pairs
will be generated. When byProbe is specified, scatter plots for one probes with nearby
20 gene pairs will be generated. When byTF is specified, scatter plot for TF expression 
and average DNA methylation at certain motif sites will be generated.
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}

\examples{
load(system.file("extdata","mee.example.rda",package = "ELMER"))
scatter.plot(mee,byProbe=list(probe=c("cg19403323"),geneNum=20),
           category="TN", save=FALSE)
scatter.plot(mee,byProbe=list(probe=c("cg19403323"),geneNum=20),
           category="TN", save=TRUE) ## save to pdf
# b. generate one probe-gene pair
scatter.plot(mee,byPair=list(probe=c("cg19403323"),gene=c("ID255928")),
            category="TN", save=FALSE,lm_line=TRUE)
}

