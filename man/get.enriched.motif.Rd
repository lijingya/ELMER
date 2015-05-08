\name{get.enriched.motif}
\alias{get.enriched.motif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
get.enriched.motif to identify the overrepresented motifs in a set of probes (HM450K)
regions.
}
\description{
get.enriched.motif is a function make use of Probes.motif data from \pkg{ELMER.data}  
package to calculate the motif enrichment Odds Ratio and  95\% confidence interval for
a given set of probes. If save is TURE, two output files will be saved: 
getMotif.XX.enriched.motifs.rda and getMotif.XX.motif.enrichment.csv (see detail).
}
\usage{
get.enriched.motif(probes.motif, probes, background.probes, lower.OR = 1.1, min.incidence = 10, 
dir.out = "./", label = NULL, save=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes.motif}{
%%     ~~Describe \code{probes.motif} here~~
A matrix contains motifs occurrence within probes regions. Probes.motif in \pkg{ELMER.data}
will be used if probes.motif is missing (detail see \code{\link{Probes.motif}}).
}
  \item{probes}{
%%     ~~Describe \code{probes} here~~
A vector lists the name of probes to define the set of probes in which motif enrichment
OR and confidence interval will be calculated.
}
  \item{background.probes}{
%%     ~~Describe \code{background.probes} here~~
A vector lists name of probes which are considered as background for motif.enrichment 
calculation (see detail).
}
  \item{lower.OR}{
A number specifies the smallest lower boundary of 95\% confidence interval for Odds Ratio.
The motif with higher lower boudnary of 95\% confidence interval for Odds Ratio than the number 
are the significantly enriched motifs (detail see reference).
}
  \item{min.incidence}{
A non-negative integer specifies the minimum incidence of motif in the given probes set. 
10 is default.
}
  \item{dir.out}{
A path specifies the directory for outputs. Default is current directory
}
  \item{label}{
A character labels the outputs such as "hypo", "hyper"
}
  \item{save}{
If save is TURE, two files will be saved: getMotif.XX.enriched.motifs.rda and 
getMotif.XX.motif.enrichment.csv (see detail).
}
}
\details{
  background.probes:
   For enhancer study, it is better to use probes within distal enhancer probes as
   background.probes. For promoter study, it is better to use probes within promoter
   regions as background.probes. Because enhancer and promoter have different CG content
   and harbors different clusters of TFs motif.
  
  save:
   if save is TRUE, two files will be save on the disk. The first file is 
   getMotif.XX.motif.enrichment.csv (XX depends on option label). This file reports 
   the Odds Ratio and 95\% confidence interval for these Odds Ratios which pass the 
   signficant cutoff (lower.OR and min.incidence). The second file is 
   getMotif.XX.enriched.motifs.rda (XX depends on option lable). This file contains
   a list R object with enriched motifs as name and probes containing the enriched 
   motif as contents. This object will be used in \code{\link{get.TFs}} function.
   if save is FALSE, the function will return a R object which is the same with second file.
  
}
\value{
A list (R object) with enriched motifs as name and probes containing the enriched 
motif as contents. And hypo.motif.enrichment.pdf plot will be generated.
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
probes <- c("cg00329272","cg10097755","cg08928189", "cg17153775","cg21156590",
"cg19749688","cg12590404","cg24517858","cg00329272","cg09010107",
"cg15386853", "cg10097755", "cg09247779","cg09181054","cg19371916")
load(system.file("extdata","mee.example.rda",package = "ELMER"))
bg <- rownames(getMeth(mee))
enriched.motif <- get.enriched.motif(probes=probes,background.probes = bg,
min.incidence=2, label="hypo")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
