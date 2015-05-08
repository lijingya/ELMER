
\name{getTCGA}
\alias{getTCGA}
\title{getTCGA to download DNA methylation, RNA expression and clinic data for 
all samples of certain cancer type from TCGA.}
\usage{
getTCGA(disease, Meth = TRUE, RNA = TRUE, Clinic = TRUE,
  basedir = "./Data", RNAtype = "gene", Methfilter = 0.2)
}
\arguments{
\item{disease}{A character specifies the disease to download in TCGA such as BLCA}

\item{Meth}{A logic if TRUE HM450K DNA methylation data will download.}

\item{RNA}{A logic if TRUE RNA-seq Hiseq-V2 from TCGA level 3 will be download.}

\item{Clinic}{A logic if TRUE clinic data will be download for that disease.}

\item{basedir}{A path shows where the data will be stored.}

\item{RNAtype}{A charactor to specify whether use isoform level or gene level.
When RNAtype=gene, gene level gene expression will be used. When isoform,
then isoform data will be used.}

\item{Methfilter}{A number. For each probe, the percentage of NA among the
all the samples should smaller than Methfilter.}
}
\value{
Download DNA methylation (HM450K)/RNAseq(HiseqV2)/Clinic data for
the specified disease from TCGA. 
}
\description{
getTCGA is a function to download DNA methylation, RNA expression and clinic data for 
all samples of certain cancer type from TCGA website. And downloaded data will be 
transform to matrixes or data frame for further analysis.
}

\references{
\url{https://tcga-data.nci.nih.gov/tcga/tcgaHome2.jsp}
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}


\examples{
getTCGA("BRCA",Meth=FALSE, RNA=FALSE, Clinic=TRUE, basedir="~")
}

