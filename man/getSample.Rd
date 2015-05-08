
\docType{methods}
\name{getSample}
\alias{getSample}
\alias{getSample,MEE.data-method}
\title{getSample to extract sample information from MEE.data object.}
\description{
getSample is a function to easily extract sample information from MEE.data object.
By specifying ID for samples, only that set of samples' informtion will be extracted.
When certain columns of data need to be extracted, just specify columns names in cols
option and sample information for wanted columns will be reported.
}
\usage{
getSample(object, ID, cols)

\S4method{getSample}{MEE.data}(object, ID, cols)
}
\arguments{
\item{object}{MEE.data object}

\item{ID}{A vector of sample ID. When specified, sample informtion only for
these samples will be output.}

\item{cols}{A vector of column names of sampleInfo slots of MEE.data object.}
}
\value{
Sample information.
}
\author{
Lijing Yao (maintainer: lijingya@usc.edu)
}


\examples{
SampleInfo <- data.frame(ID=paste0("sample",1:5),
TN=c("Tumor","Tumor","Normal","Normal","Tumor"))
mee <- fetch.mee(sample = SampleInfo)
Samples <- getSample(mee,ID = "sample2") ## get sample2's information
Samples <- getSample(mee, cols = "TN")  ## get 'TN' information for each samples
}

