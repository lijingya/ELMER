
\name{getClinic}
\alias{getClinic}
\title{getClinic to download clinic data for certain cancer types from TCGA website.}

\description{
getClinic is a function to download latest version of clinic data for 
all samples of certain cancer types from TCGA website.
}

\usage{
getClinic(disease, basedir = "./Data")
}
\arguments{
\item{disease}{A character specifies the disease to download from TCGA such as BLCA}

\item{basedir}{A path shows where the data will be stored.}
}
\value{
Download all clinic information for the specified disease.
}

\references{
\url{https://tcga-data.nci.nih.gov/tcga/tcgaHome2.jsp}
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}


