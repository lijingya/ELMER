## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
BiocStyle::latex()

## ----installing, eval=FALSE-------------------------------------------------------------
#  install.packages(devtools)
#  library(devtools);
#  devtools::install_github("lijingya/ELMER");

## ----example.data.run, echo=FALSE-------------------------------------------------------
if(!file_test("-d", "./ELMER.example")) {
  URL <- "https://dl.dropboxusercontent.com/u/61961845/ELMER.example.tar.gz"
  download.file(URL,destfile = "ELMER.example",method= "wget",
                extra = c("--no-check-certificate -a download.log"))
  untar("./ELMER.example")
}
library(ELMER)

## ----example.data, eval=FALSE-----------------------------------------------------------
#  #Example file download from URL: https://dl.dropboxusercontent.com/u/61961845/ELMER.example.tar.gz
#  URL <- "https://dl.dropboxusercontent.com/u/61961845/ELMER.example.tar.gz"
#  download.file(URL,destfile = "ELMER.example",method= "wget",
#                extra = c("--no-check-certificate -a download.log"))
#  untar("./ELMER.example")
#  library(ELMER)

## ----tcga.pipe, cache=TRUE--------------------------------------------------------------
TCGA.pipe("LUSC",wd="./ELMER.example",cores=detectCores()/2,permu.size=300,
          analysis = c("distal.enhancer","diffMeth","pair","motif","TF.search"),
          diff.dir="hypo",rm.chr=paste0("chr",c(1:22,"X","Y")))

## ----dna.methylation.data, cache=TRUE---------------------------------------------------
load("./ELMER.example/Result/LUSC/LUSC_meth_refined.rda")
Meth[1:10, 1:2]

## ----gene.expression.data, cache=TRUE---------------------------------------------------
load("./ELMER.example/Result/LUSC/LUSC_RNA_refined.rda")
GeneExp[1:10, 1:2]

## ----sample.information, cache=TRUE-----------------------------------------------------
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=T)
head(getSample(mee))

## ----probe.information, cache=TRUE------------------------------------------------------
probe <- ReadBed(system.file("extdata","Illumina-methyl-450K-manifest.hg19.bed.xz",
                             package = "ELMER"))
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=T, probeInfo=probe)
getProbeInfo(mee)

## ----gene.information, cache=TRUE-------------------------------------------------------
load(system.file("extdata","UCSC_gene_hg19.rda",package = "ELMER"))
## In TCGA expression data, geneIDs were used as the rowname for each row. However, numbers 
## can't be the rownames, "ID" was added to each gene id functioning as the rowname.
## If your geneID is consistent with the rownames of the gene expression matrix, adding "ID" 
## to each geneID can be skipped.
txs$GENEID <- paste0("ID",txs$GENEID)            
geneInfo <- promoters(txs,upstream = 0, downstream = 0)
save(geneInfo,file="./ELMER.example/Result/LUSC/geneAnnot.rda")
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=T, geneInfo=txs)
getGeneInfo(mee)

## ----mee.data, cache=TRUE---------------------------------------------------------------
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=T, probeInfo=probe, geneInfo=txs)
mee

## ----selection.of.probes.within.biofeatures, cache=TRUE---------------------------------
#get distal enhancer probes that are 2kb away from TSS and overlap with REMC and FANTOM5 
#enhancers on chromosome 1
Probe <- get.feature.probe(probe=probe, rm.chr=paste0("chr",c(2:22,"X","Y")))
save(Probe,file="./ELMER.example/Result/LUSC/probeInfo_feature.rda")

## ----identifying.differentially.methylated.probes, cache=TRUE---------------------------
## fetch.mee can take path as input.
mee <- fetch.mee(meth="./ELMER.example/Result/LUSC/LUSC_meth_refined.rda",
                 exp="./ELMER.example/Result/LUSC/LUSC_RNA_refined.rda", TCGA=T, 
                 probeInfo="./ELMER.example/Result/LUSC/probeInfo_feature.rda", 
                 geneInfo="./ELMER.example/Result/LUSC/geneAnnot.rda") 

sig.diff <- get.diff.meth(mee, cores=detectCores()/2, dir.out ="./ELMER.example/Result/LUSC", 
                          diff.dir="hypo", pvalue = 0.01)

sig.diff$hypo[1:10,]   ## significantly hypomethylated probes

# get.diff.meth automatically save output files. 
# getMethdiff.hypo.probes.csv contains statistics for all the probes.
# getMethdiff.hypo.probes.significant.csv contains only the significant probes.
dir(path = "./ELMER.example/Result/LUSC", pattern = "getMethdiff")  

## ----identifying.putative.probe.gene.pairs, cache=TRUE----------------------------------
### identify target gene for significantly hypomethylated probes.

Sig.probes <- read.csv("./ELMER.example/Result/LUSC/getMethdiff.hypo.probes.significant.csv",
                       stringsAsFactors=F)[,1]  
head(Sig.probes)  # significantly hypomethylated probes

## Collect nearby 20 gene for Sig.probes
nearGenes <-GetNearGenes(TRange=getProbeInfo(mee,probe=Sig.probes),
                         geneAnnot=getGeneInfo(mee),cores=detectCores()/2)

## Identify significant probe-gene pairs
Hypo.pair <-get.pair(mee=mee,probes=Sig.probes,nearGenes=nearGenes,
                     permu.dir="./ELMER.example/Result/LUSC/permu",permu.size=300,Pe = 0.01,
                     dir.out="./ELMER.example/Result/LUSC",cores=detectCores()/2,label= "hypo")

head(Hypo.pair)  ## significant probe-gene pairs

# get.pair automatically save output files. 
#getPair.hypo.all.pairs.statistic.csv contains statistics for all the probe-gene pairs.
#getPair.hypo.pairs.significant.csv contains only the significant probes.
dir(path = "./ELMER.example/Result/LUSC", pattern = "getPair") 

## ----motif.enrichment.analysis.on.selected.probes, cache=TRUE---------------------------
### identify enriched motif for significantly hypomethylated probes which 
##have putative target genes.

Sig.probes.paired <- read.csv("./ELMER.example/Result/LUSC/getPair.hypo.pairs.significant.csv",
                              stringsAsFactors=F)[,1]  
head(Sig.probes.paired) # significantly hypomethylated probes with putative target genes

enriched.motif <-get.enriched.motif(probes=Sig.probes.paired, 
                                    dir.out="./ELMER.example/Result/LUSC", label="hypo",
                                    min.incidence = 10,lower.OR = 1.1)
names(enriched.motif)  # enriched motifs

# get.enriched.motif automatically save output files. 
# getMotif.hypo.enriched.motifs.rda contains enriched motifs and the probes with the motif. 
# getMotif.hypo.motif.enrichment.csv contains summary of enriched motifs.
dir(path = "./ELMER.example/Result/LUSC", pattern = "getMotif") 

# motif enrichment figure will be automatically generated.
dir(path = "./ELMER.example/Result/LUSC", pattern = "motif.enrichment.pdf") 

## ----identifying.regulatory.tf, cache=TRUE----------------------------------------------
### identify regulatory TF for the enriched motifs

load("./ELMER.example/Result/LUSC/getMotif.hypo.enriched.motifs.rda")
TF <- get.TFs(mee=mee, enriched.motif=enriched.motif,dir.out="./ELMER.example/Result/LUSC", 
              cores=detectCores()/2, label= "hypo")

# get.TFs automatically save output files. 
# getTF.hypo.TFs.with.motif.pvalue.rda contains statistics for all TF with average 
# DNA methylation at sites with the enriched motif.
# getTF.hypo.significant.TFs.with.motif.summary.csv contains only the significant probes.
dir(path = "./ELMER.example/Result/LUSC", pattern = "getTF")  

# TF ranking plot based on statistics will be automatically generated.
dir(path = "./ELMER.example/Result/LUSC/TFrankPlot", pattern = "pdf") 

## ----figure1, fig.cap="\\label{fig:cg19403323.byprobe} Each scatter plot shows the methylation level of an example probe cg19403323 in all LUSC samples plotted against the expression of one of 20 adjacent genes.", fig.scap="Scatter Plot of 20 nearby genes", cache=TRUE, fig.height=5.5----
scatter.plot(mee,byProbe=list(probe=c("cg19403323"), geneNum=20), 
             category="TN", save=FALSE)

## ----figure2, fig.cap="\\label{fig:cg19403323_SYT14.bypair} Scatter plot shows the methylation level of an example probe cg19403323 in all LUSC samples plotted against the expression of the putative target gene SYT14.", fig.scap="Scatter Plot of One Pair", cache=TRUE, fig.height=3.5----
scatter.plot(mee, byPair=list(probe=c("cg19403323"), gene=c("ID255928")),
             category="TN", save=FALSE, lm_line=TRUE)

## ----figure3, fig.cap="\\label{fig:TP53_TP63_TP73.byTF} Each scatter plot shows the average methylation level of sites with the TP53 motif in all LUSC samples plotted against the expression of the transcription factor TP53, TP63, TP73 respectively.", fig.scap="TF expression vs. average DNA methylation", cache=TRUE, fig.height=3.5----
load("ELMER.example/Result/LUSC/getMotif.hypo.enriched.motifs.rda")
scatter.plot(mee,byTF=list(TF=c("TP53","TP63","TP73"),
                           probe=enriched.motif[["TP53"]]), category="TN",
             save=FALSE, lm_line=TRUE)

## ----makepair, cache=TRUE---------------------------------------------------------------
# Make a "Pair" object for schematic.plot
pair <- fetch.pair(pair="./ELMER.example/Result/LUSC/getPair.hypo.pairs.significant.withmotif.csv",
                   probeInfo = "./ELMER.example/Result/LUSC/probeInfo_feature.rda",
                   geneInfo = "./ELMER.example/Result/LUSC/geneAnnot.rda")

## ----figure4, fig.cap="\\label{fig:cg19403323.schematic.byProbe} The schematic plot shows probe colored in blue and the location of nearby 20 genes. The genes significantly linked to the probe were shown in red.", fig.scap="Nearby Genes", cache=TRUE, fig.height=3.5----
grid::grid.newpage()
schematic.plot(pair=pair, byProbe="cg19403323", save=FALSE)

## ----figure5, fig.cap="\\label{fig:ID255928.schematic.byGene} The schematic plot shows the gene colored in red and all blue colored probes, which are significantly linked to the expression of this gene.", fig.scap="Nearby Probes", cache=TRUE, fig.height=3.5----
grid::grid.newpage()
schematic.plot(pair=pair, byGene="ID255928", save=FALSE)

## ----figure6, fig.cap="\\label{fig:hypo.motif.enrichment} The plot shows the Odds Ratio (x axis) for the selected motifs with OR above 1.3 and lower boundary of OR above 1.3. The range shows the 95\\% confidence interval for each Odds Ratio.", fig.scap="Motif Enrichment", cache=TRUE, fig.height=2----
motif.enrichment.plot(motif.enrichment="./ELMER.example/Result/LUSC/getMotif.hypo.motif.enrichment.csv", 
                      significant=list(OR=1.3, lowerOR=1.3), label="hypo", save=FALSE)

## ----figure7, fig.cap="\\label{fig:TP53.TFrankPlot} Shown are TF ranking plots based on the score (-log(P value)) of association between TF expression and DNA methylation of the TP53 motif in the LUSC cancer type . The dashed line indicates the boundary of the top 5\\% association score. The top 3 associated TFs and the TF family members (dots in red) that are associated with that specific motif are labeled in the plot.", fig.scap="TF Ranking", cache=TRUE, fig.height=4.5----
load("./ELMER.example/Result/LUSC/getTF.hypo.TFs.with.motif.pvalue.rda")
TF.rank.plot(motif.pvalue=TF.meth.cor, motif="TP53", TF.label=list(TP53=c("TP53","TP63","TP73")),
             save=FALSE)

## ----finalsession, cache=TRUE-----------------------------------------------------------
sessionInfo()

