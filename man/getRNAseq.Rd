
\name{getRNAseq}
\alias{getRNAseq}
\title{getRNAseq to download all RNAseq data for a certain cancer type from TCGA.}
\description{
getRNAseq is a function to download RNAseq data for all samples of a certain cancer
type from TCGA
}
\usage{
getRNAseq(disease, basedir = "./Data")
}
\arguments{
\item{disease}{A character specifies disease in TCGA such as BLCA}

\item{basedir}{A path shows where the data will be stored.}
}
\value{
Download all RNA seq level 3 data for the specified disease.
}


\references{
\url{https://tcga-data.nci.nih.gov/tcga/tcgaHome2.jsp}
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}
