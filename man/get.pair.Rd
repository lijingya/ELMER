\name{get.pair}
\alias{get.pair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.pair to predict enhancer-gene linkages.
}
\description{
get.pair is a function to predict enhancer-gene linkages using associations between
DNA methylation at enhancer CpG sites and expression of 20 nearby genes of the CpG sites
(see reference). Two files will be saved if save is true: getPair.XX.all.pairs.statistic.csv
and getPair.XX.pairs.significant.csv (see detail).
}
\usage{
get.pair(mee, probes, nearGenes, percentage = 0.2, permu.size = 10000, permu.dir = NULL, 
Pe = 0.001, dir.out = "./", diffExp = FALSE, cores = NULL, portion=0.3,   
label = NULL, save=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mee}{
A MEE.data object contains at least meth, exp, probeInfo, geneInfo slots.
}
  \item{probes}{
A vector lists name of probes that need to be linked to genes.
}
  \item{nearGenes}{
A list (R object) containing output of GetNearGenes function or a path of XX.rda file
containing output of GetNearGenes function.
}
  \item{percentage}{
A number ranges from 0 to 1 specifying the percentage of samples of methylated and 
unmethylated groups used to link probes to genes. Default is 0.2.
}
  \item{permu.size}{
A number specifies the number of permuation. Default is 1000.
}
  \item{permu.dir}{
A path shows the directory of permutation outputs.
}
  \item{Pe}{
A number specifies the empircal pvalue cutoff for defining signficant pairs.
Default is 0.01
}
  \item{portion}{
A number specify the cut point to define binary methlation level for probe loci.
Default is 0.3. When beta value is above 0.3, the probe is methylated and 
vice versa. For one probe, the percentage of methylated or unmethylated samples 
should be above 0.05.
}

  \item{dir.out}{
A path specifies the directory for outputs of get.pair function. Default is current directory
}
  \item{diffExp}{
A logic. Default is FALSE. If TRUE, t test will be applied to test whether putative 
target gene are differentially expressed between two control and experimental 
groups.
}
  \item{cores}{
A interger which defines the number of cores to be used in parallel process.
Default is NULL: no parallel process.
}
  \item{label}{
A character labels the outputs.
}
  \item{save}{
A logic. If save is true, two files will be saved for publication or analysis re-usage 
purpose: getPair.XX.all.pairs.statistic.csv and getPair.XX.pairs.significant.csv (see detail)
}
}
\details{
save:
 When save is TRUE, function will generate two XX.csv files.
 The first one is named getPair.XX.all.pairs.statistic.csv ( XX depends on option label). 
 This file contains all statistic results for each probe-gene pair. Based on this file, user 
 can change different P value or sig.dir cutoff to select the significant results
 without redo the analysis. The second file is named getPair.XX.pairs.significant.csv
 (XX depends on option label). This file contains statistic results for the probes 
 that pass the significant criteria (Pe).
 When save is FALSE, a data frame R object will be generate which contains the same
 information with the second file.

}
\value{
A data frame contains statistic result for significant pairs
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
nearGenes <-GetNearGenes(TRange=getProbeInfo(mee,probe=c("cg00329272","cg10097755")),
                         geneAnnot=getGeneInfo(mee))
                         Hypo.pair <-get.pair(mee=mee,probes=c("cg00329272","cg10097755"),
                                              nearGenes=nearGenes,permu.size=5,Pe = 0.2,
                                              dir.out="./",
                                              label= "hypo")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
