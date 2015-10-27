\name{get.permu}
\alias{get.permu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.permu to generate permutation results for calculation of empirical P values for
each enhancer-gene linkage.
}
\description{
get.permu is a function to use the same statistic model to calculate random enhancer-gene
pairs. Based on the permutation value, empirical P value can be calculated for the
real enhancer-gene pair (see reference).
}
\usage{
get.permu(mee, geneID, percentage = 0.2, rm.probes = NULL,
portion = 0.3, permu.size = 10000, permu.dir = NULL, cores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mee}{
A MEE.data object contains at least meth, exp, probeInfo, geneInfo.
}
  \item{geneID}{
A vector lists gene id which need to be have permutation.
}
  \item{percentage}{
A number ranges from 0 to 1 specifying the percentage of samples of control and experimental
groups used to link probes to genes. Default is 0.2.
}
  \item{portion}{
A number specify the cut point to define binary methlation level for probe loci.
Default is 0.3. When beta value is above 0.3, the probe is methylated and 
vice versa. For one probe, the percentage of methylated or unmethylated samples 
should be above 0.05.
}

  \item{rm.probes}{
A vector lists name of probes which belongs to the set of probes fed into \code{\link{get.pair}}
function.
}
  \item{permu.size}{
A number specifies the number of permuation. Default is 1000.
}
  \item{permu.dir}{
A path shows the directory of permutation outputs
}
  \item{cores}{
A interger which defines the number of cores to be used in parallel process.
Default is NULL: no parallel process.
}
}

\value{
Certain number of permutation for each gene of interets.
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
\note{
Permutation is the most time consuming step. It is recommended to use multiple 
cores for this step. Default permutation time is 1000 which may need 12 hrs by 4 cores. 
However 10,000 permutations is recommended to get high confidence results. But it 
may cost 2 days.
}


\examples{
load(system.file("extdata","mee.example.rda",package = "ELMER"))
permu <-get.permu(mee=mee,geneID=rownames(getExp(mee)),
                  rm.probes=c("cg00329272","cg10097755"),
                  permu.size=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
