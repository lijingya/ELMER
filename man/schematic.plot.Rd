
\name{schematic.plot}
\alias{schematic.plot}
\title{schematic.plot to plot schematic plots showing the locations of genes and probes.}
\usage{
schematic.plot(pair, byProbe, byGene, byCoordinate = list(chr = c(), start =
  c(), end = c()), dir.out = "./", save = TRUE, ...)
}
\arguments{
\item{pair}{A Pair object. All slots of Pair class should be included}

\item{byProbe}{A vector of probe names.}

\item{byGene}{A vector of gene ID}

\item{byCoordinate}{A list contains chr, start and end.
byCoordinate=list(chr=c(),start=c(),end=c()).}

\item{dir.out}{A path specifies the directory for outputs. Default is current directory}

\item{save}{A logic. If true, figures will be saved to dir.out.}

\item{...}{Parameters for GetNearGenes function. See \code{\link{GetNearGenes}}}
}
\value{
Schematic plots will be produced.
}
\description{
schematic.plot is a function to plot schematic plots showing the locations of genes and probes.
}
\details{
byProbes:
When a vector of probes' name are provided,
function will produce schematic plots for each individual probes.
The schematic plot contains probe, nearby 20 (or the number of gene user specified.)
genes and the significantly linked gene to the probe.
  
byGene:
When a vector of gene ID are provided, function will produce schematic plots 
for each individual genes. The schematic plot contains the gene and all the 
significantly linked probes.

byCoordinate:
When a genomic coordinate is provided, function will 
produce a schematic plot for this coordinate. The schematic plot contains 
all genes and significantly linked probes in the range and the significant links.
  
 
}

\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}

\examples{
library(grid)
load(system.file("extdata","mee.example.rda",package = "ELMER"))
nearGenes <-GetNearGenes(TRange=getProbeInfo(mee,probe=c("cg00329272","cg19403323")),
                        geneAnnot=getGeneInfo(mee))
Hypo.pair <-get.pair(mee=mee,probes=c("cg00329272","cg19403323"),
                    nearGenes=nearGenes,permu.size=5,Pe = 0.2,dir.out="./",
                    label= "hypo")
pair <- fetch.pair(pair=Hypo.pair,
                  probeInfo = getProbeInfo(mee),
                  geneInfo = getGeneInfo(mee))
# a. generate schematic plot of one probe with nearby 20 genes and label
#the gene significantly linked with the probe.
grid.newpage()
schematic.plot(pair=pair, byProbe="cg19403323" ,save=FALSE)
#b. generate schematic plot of ont gene with the probe which the gene significanlty linked to.
grid.newpage()
schematic.plot(pair=pair, byGene="ID255928",save=FALSE)
}

