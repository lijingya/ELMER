\name{get450K}
\alias{get450K}
\title{get450K to download HM40K DNA methylation data for certain cancer types 
from TCGA website.}

\description{
get450K is a function to download latest version of HM450K DNA methylation for 
all samples of certain cancer types from TCGA website.
}

\usage{
get450K(disease, basedir = "./Data")
}
\arguments{
\item{disease}{A character specifies the disease to download from TCGA such as BLCA}

\item{basedir}{A path shows where the data will be stored.}
}
\value{
Download all DNA methylation from HM450K level 3 data for the
specified disease.
}

\references{
\url{https://tcga-data.nci.nih.gov/tcga/tcgaHome2.jsp}
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}


