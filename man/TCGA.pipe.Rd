
\name{TCGA.pipe}
\alias{TCGA.pipe}
\title{ELMER analysis pipeline for TCGA data.}
\usage{
TCGA.pipe(disease, analysis = "all", wd = "./", cores = NULL,
  Data = NULL, diff.dir = "hypo", ...)
}
\arguments{
\item{disease}{TCGA disease name in short form such as COAD}

\item{analysis}{A vector of characters listing the analysis need to be done.
Analysis can be "download","distal.probes","diffMeth","pair","motif","TF.search".
Default is "all" meaning all the analysis will be processed.}

\item{wd}{A path shows working dirctory. Default is "./"}

\item{cores}{A interger defines number of core to be used in parallel process.
Default is NULL: don't use parallel process.}

\item{Data}{A path shows the folder containing DNA methylation, expression and clinic data}

\item{diff.dir}{A character can be "hypo" or "hyper", showing dirction DNA methylation changes.  
If it is "hypo", get.diff.meth function will identify all significantly hypomethylated
CpG sites; If "hyper", get.diff.meth function will identify all significantly hypermethylated
CpG sites}

\item{...}{A list of parameters for functions: GetNearGenes, get.feature.probe,
get.diff.meth, get.pair,}
}
\value{
Different analysis results.
}
\description{
ELMER analysis pipeline for TCGA data. This pipeline combine every steps of \pkg{ELMER}
analyses: get.feature.probe, get.diff.meth, get.pair, get.permu, get.enriched.motif and get.TFs.
Every steps' results are saved.
}
\examples{
\dontrun{
distal.probe <- TCGA.pipe(disease = "LUSC", analysis="Probe.selection", wd="~/")
}
}

