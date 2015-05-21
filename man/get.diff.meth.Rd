\name{get.diff.meth}
\alias{get.diff.meth}
\title{
get.diff.meth to identify hypo/hyper-methylated CpG sites on HM450K between 
control and experimental groups such as normal verus tumor samples.
}
\description{
get.diff.meth applys one-way t-test to identify the CpG sites that are significantly 
hypo/hyper-methyalated using proportional samples (defined by percentage option) 
from control and experimental groups. The P values will be adjusted by Benjamini-Hochberg 
method. Option pvalue and sig.dif will be the criteria (cutoff) for selecting significant 
differentially methylated CpG sites. If save is TURE, two getMethdiff.XX.csv files
will be generated (see detail).
}
\usage{
get.diff.meth(mee, diff.dir = "hypo", cores = NULL, percentage = 0.2, pvalue = 0.01, sig.dif = 0.3, dir.out = "./", save = TRUE)
}
\arguments{
  \item{mee}{
 A MEE.data class object contains at least methy and probeInfo slots.
}
  \item{diff.dir}{
A character can be "hypo" or "hyper", showing dirction DNA methylation changes.  
If it is "hypo", get.diff.meth function will identify all significantly hypomethylated
CpG sites; If "hyper", get.diff.meth function will identify all significantly hypoermethylated
CpG sites
}
  \item{cores}{
A interger which defines the number of cores to be used in parallel process.
Default is NULL: no parallel process.
}
  \item{percentage}{
A number ranges from 0 to 1 specifying the percentage of samples from control and 
experimental groups that are used to identify the differential methylation. 
Default is 0.2.
}
  \item{pvalue}{
A number specifies the significant P value (adjusted P value by BH) cutoff for 
selecting significant hypo/hyper-methylated probes. Default is 0.01.
}
  \item{sig.dif}{
A number specifies the smallest DNA methylation difference as a cutoff for selecting 
significant hypo/hyper-methylated probes. Default is 0.3.
}
  \item{dir.out}{
A path specifies the directory for outputs. Default is is current directory.
}
  \item{save}{
A logic. When TRUE, two getMethdiff.XX.csv files will be generated (see detail)
}
}
\details{
save:
 When save is TRUE, function will generate two XX.csv files.The first one is 
 named getMethdiff.hypo.probes.csv (or getMethdiff.hyper.probes.csv depends
 on diff.dir). The first file contains all statistic results for each probe. Based on this
 file, user can change different P value or sig.dir cutoff to select the significant results
 without redo the analysis. The second file is named getMethdiff.hypo.probes.significant.csv
 (or getMethdiff.hyper.probes.significant.csv depends on diff.dir). This file contains
 statistic results for the probes that pass the significant criteria (P value and sig.dir).
 When save is FALSE, a data frame R object will be generate which contains the same
 information with the second file.

}
\value{
A data frame contains statistics from differential analysis for each probes.
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
Hypo.probe <- get.diff.meth(mee, diff.dir="hypo") # get hypomethylated probes
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
