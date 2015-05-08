## ELMER
### An R/Bioconductor Tool Inferring Regulatory Element Landscapes and Transcription Factor Networks Using Methylomes

#### Installing and loading ELMER
To obtain a copy of ELMER, you will need to install devtools and ELMER.data which contains essential data for running ELMER package 


```r
install.packages(devtools)
library(devtools);
devtools::install_github("lijingya/ELMER.data");
devtools::install_github("lijingya/ELMER");
```
Then you can load the package and see an introduction with
```r
library(ELMER)
openVignette("ELMER")
```
```
# Please select a vignette:
# 1: ELMER - ELMER: Inferring Regulatory Element Landscapes and Transcription Factor Networks Using Methylomes
# Selection: 1
```
Or you can have the vignette, and sourcecode for the vignette open directly in the browser
```r
browseVignettes("ELMER")
```
