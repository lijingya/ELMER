\name{get.TFs}
\alias{get.TFs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.TFs to identify regulatory TFs.
}
\description{
get.TFs is a function to identify regulatory TFs based on motif analysis and association
analysis between the probes containing a particular motif and expression of all known TFs.
If save is ture, two files will be saved: getTF.XX.significant.TFs.with.motif.summary.csv
and getTF.hypo.TFs.with.motif.pvalue.rda (see detail).
}
\usage{
get.TFs(mee, enriched.motif, TFs, motif.relavent.TFs, percentage = 0.2, 
dir.out = "./", label = NULL, cores = NULL,save=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mee}{
A MEE.data object contains at least meth, exp, probeInfo, geneInfo.
}
  \item{enriched.motif}{
A list containing output of \code{\link{get.enriched.motif}} function or a 
path of XX.rda file containing output of get.enriched.motif function.
}
  \item{TFs}{
A data.frame containing TF GeneID and Symbol or a path of XX.csv file containing 
TF GeneID and Symbol. If missing, human.TF list will be used (human.TF data in
\pkg{ELMER.data}). For detail information, refer the reference paper.
}
  \item{motif.relavent.TFs}{
A list containing motif as names and relavent TFs as contents for each list element or a 
path of XX.rda file containing a list as above. If missing, motif.relavent.TFs will
be used (motif.relavent.TFs data in \pkg{ELMER.data}). For detail information, 
refer the reference paper.
}
  \item{percentage}{
A number ranges from 0 to 1 specifying the percentage of samples of control and experimental
groups used to link probes to genes. Default is 0.2.
}
  \item{dir.out}{
A path specifies the directory for outputs of get.pair function. Default is current directory
}
  \item{cores}{
A interger which defines the number of cores to be used in parallel process.
Default is NULL: no parallel process.
}
  \item{label}{
A character labels the outputs.
}
  \item{save}{
A logic. If save is ture, two files will be saved: getTF.XX.significant.TFs.with.motif.summary.csv
and getTF.hypo.TFs.with.motif.pvalue.rda (see detail). If save is false, a data frame
contains the same content with the first file.
}
}
\details{
save:
  If save is ture, two files will be saved. The first file is
  getTF.XX.significant.TFs.with.motif.summary.csv (XX depends on option lable). 
  This file contain the regulatory TF significantly associate with average DNA methylation
  at particular motif sites. The second file is getTF.hypo.TFs.with.motif.pvalue.rda
  (XX depends on option label). This file contains a matrix storing the statistic 
  results for significant associations between TFs (row) and average DNA methylation 
  at motifs (column). If save is false, a data frame which contains the same content 
  with the first file will be reported.
  
}
\value{
Potential responsible TFs will be reported.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
load(system.file("extdata","mee.example.rda",package = "ELMER"))
enriched.motif <- list("TP53"= c("cg00329272", "cg10097755", "cg08928189",
                                 "cg17153775", "cg21156590", "cg19749688", "cg12590404",
                                 "cg24517858", "cg00329272", "cg09010107", "cg15386853",
                                 "cg10097755", "cg09247779", "cg09181054"))
TF <- get.TFs(mee, enriched.motif,
              TFs=data.frame(GeneID=c("ID7157","ID8626","ID7161"),
              Symbol=c("TP53","TP63","TP73"),
              stringsAsFactors = FALSE),
              label="hypo")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
