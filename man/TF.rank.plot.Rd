
\name{TF.rank.plot}
\alias{TF.rank.plot}
\title{TF.rank.plot to plot the scores (-log10(P value)) which assess the
correlation between TF expression and average DNA methylation at motif sites.}
\usage{
TF.rank.plot(motif.pvalue, motif, TF.label, dir.out = "./", save = TRUE)
}
\arguments{
\item{motif.pvalue}{A matrix or a path specifying location of "XXX.with.motif.pvalue.rda"
which is output of get.TFs.}

\item{motif}{A vector of charactors specify the motif to plot}

\item{TF.label}{A list shows the labels for each motif. If TF.label is not specified,
the motif relevant TF and top3 TF will be labeled.}

\item{dir.out}{A path specify the directory to which the figures will be saved.
Current directory is default.}

\item{save}{A logic. If true (default), figure will be saved to dir.out.}
}
\value{
A plot shows the score (-log(P value)) of association between TF
expression and DNA methylation at sites of a certain motif.
}
\description{
TF.rank.plot is a function to plot the scores (-log10(P value)) which assess the
correlation between TF expression and average DNA methylation at motif sites. The the motif 
relevant TF and top3 TFs will be labeled in a different color.
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}
\examples{
load(system.file("extdata","getTF.hypo.TFs.with.motif.pvalue.rda",package="ELMER"))
TF.rank.plot(motif.pvalue=TF.meth.cor, motif="TP53", TF.label=list(TP53=c("TP53","TP63","TP73")),
            save=FALSE)
}

