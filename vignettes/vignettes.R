## ----installing, eval=FALSE,hide=TRUE, message=FALSE,include=FALSE-------
install.packages(devtools)
library(devtools);
devtools::install_github("lijingya/ELMER");

## ----install, eval=FALSE-------------------------------------------------
source("http://bioconductor.org/biocLite.R")
biocLite("ELMER")

## ----example.data, eval=FALSE--------------------------------------------
## #Example file download from URL: https://dl.dropboxusercontent.com/u/61961845/ELMER.example.tar.gz
URL <- "https://dl.dropboxusercontent.com/u/61961845/ELMER.example.tar.gz"
download.file(URL,destfile = "ELMER.example.tar.gz",method= "curl")
untar("./ELMER.example.tar.gz")
library(ELMER)

## ----tcga.pipe, cache=TRUE-----------------------------------------------
TCGA.pipe("LUSC",wd="./ELMER.example",cores=detectCores()/2,permu.size=300,
          analysis = c("distal.probes","diffMeth","pair","motif","TF.search"),
          diff.dir="hypo",rm.chr=paste0("chr",c("X","Y")))

## ----dna.methylation.data, cache=TRUE------------------------------------
load("./ELMER.example/Result/LUSC/LUSC_meth_refined.rda")
Meth[1:10, 1:2]

## ----gene.expression.data, cache=TRUE------------------------------------
load("./ELMER.example/Result/LUSC/LUSC_RNA_refined.rda")
GeneExp[1:10, 1:2]

## ----sample.information, cache=TRUE--------------------------------------
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=TRUE)
head(getSample(mee))

## ----probe.information, cache=TRUE---------------------------------------
probe <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, what="Locations")
probe <- GRanges(seqnames=probe$chr,
                 ranges=IRanges(probe$pos,
                                width=1,
                                names=rownames(probe)),
                 strand=probe$strand,
                 name=rownames(probe))
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=TRUE, probeInfo=probe)
getProbeInfo(mee)

## ----gene.information, cache=TRUE----------------------------------------
geneAnnot <- txs()
## In TCGA expression data, geneIDs were used as the rowname for each row. However, numbers
## can't be the rownames, "ID" was added to each gene id functioning as the rowname.
## If your geneID is consistent with the rownames of the gene expression matrix, adding "ID"
## to each geneID can be skipped.
geneAnnot$GENEID <- paste0("ID",geneAnnot$GENEID)
geneInfo <- promoters(geneAnnot,upstream = 0, downstream = 0)
save(geneInfo,file="./ELMER.example/Result/LUSC/geneAnnot.rda")
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=TRUE, geneInfo=geneInfo)
getGeneInfo(mee)

## ----mee.data, cache=TRUE------------------------------------------------
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=TRUE, probeInfo=probe, geneInfo=geneInfo)
mee

## ----selection.of.probes.within.biofeatures, cache=TRUE------------------
#get distal enhancer probes that are 2kb away from TSS and overlap with REMC and FANTOM5
#enhancers on chromosome 1
Probe <- get.feature.probe(rm.chr=paste0("chr",c(2:22,"X","Y")))
save(Probe,file="./ELMER.example/Result/LUSC/probeInfo_feature_distal.rda")

## ----identifying.differentially.methylated.probes, cache=TRUE------------
## fetch.mee can take path as input.
mee <- fetch.mee(meth="./ELMER.example/Result/LUSC/LUSC_meth_refined.rda",
                 exp="./ELMER.example/Result/LUSC/LUSC_RNA_refined.rda", TCGA=TRUE,
                 probeInfo="./ELMER.example/Result/LUSC/probeInfo_feature_distal.rda",
                 geneInfo="./ELMER.example/Result/LUSC/geneAnnot.rda")

sig.diff <- get.diff.meth(mee, cores=detectCores()/2, dir.out ="./ELMER.example/Result/LUSC",
                          diff.dir="hypo", pvalue = 0.01)


sig.diff[1:10,]   ## significantly hypomethylated probes

# get.diff.meth automatically save output files.
# getMethdiff.hypo.probes.csv contains statistics for all the probes.
# getMethdiff.hypo.probes.significant.csv contains only the significant probes which
# is the same with sig.diff
dir(path = "./ELMER.example/Result/LUSC", pattern = "getMethdiff")

## ----identifying.putative.probe.gene.pairs, cache=TRUE-------------------
### identify target gene for significantly hypomethylated probes.

Sig.probes <- read.csv("./ELMER.example/Result/LUSC/getMethdiff.hypo.probes.significant.csv",
                       stringsAsFactors=FALSE)[,1]
head(Sig.probes)  # significantly hypomethylated probes

## Collect nearby 20 genes for Sig.probes
nearGenes <-GetNearGenes(TRange=getProbeInfo(mee,probe=Sig.probes),
                         geneAnnot=getGeneInfo(mee),cores=detectCores()/2)

## Identify significant probe-gene pairs
Hypo.pair <-get.pair(mee=mee,probes=Sig.probes,nearGenes=nearGenes,
                     permu.dir="./ELMER.example/Result/LUSC/permu",permu.size=300,Pe = 0.01,
                     dir.out="./ELMER.example/Result/LUSC",cores=detectCores()/2,label= "hypo")

head(Hypo.pair)  ## significant probe-gene pairs

# get.pair automatically save output files.
#getPair.hypo.all.pairs.statistic.csv contains statistics for all the probe-gene pairs.
#getPair.hypo.pairs.significant.csv contains only the significant probes which is
# same with Hypo.pair.
dir(path = "./ELMER.example/Result/LUSC", pattern = "getPair")

## ----motif.enrichment.analysis.on.selected.probes, cache=TRUE------------
### identify enriched motif for significantly hypomethylated probes which
##have putative target genes.

Sig.probes.paired <- read.csv("./ELMER.example/Result/LUSC/getPair.hypo.pairs.significant.csv",
                              stringsAsFactors=FALSE)[,1]
head(Sig.probes.paired) # significantly hypomethylated probes with putative target genes

enriched.motif <-get.enriched.motif(probes=Sig.probes.paired,
                                    dir.out="./ELMER.example/Result/LUSC", label="hypo",
                                    min.incidence = 10,lower.OR = 1.1)
names(enriched.motif)  # enriched motifs
head(enriched.motif["TP53"]) ## probes in the given set that have TP53 motif.

# get.enriched.motif automatically save output files.
# getMotif.hypo.enriched.motifs.rda contains enriched motifs and the probes with the motif.
# getMotif.hypo.motif.enrichment.csv contains summary of enriched motifs.
dir(path = "./ELMER.example/Result/LUSC", pattern = "getMotif")

# motif enrichment figure will be automatically generated.
dir(path = "./ELMER.example/Result/LUSC", pattern = "motif.enrichment.pdf")

## ----identifying.regulatory.tf, cache=TRUE-------------------------------
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

## ----figure1,fig.width=11, fig.height=10,fig.scap="Scatter Plot of 20 nearby genes",cache=TRUE----
scatter.plot(mee,byProbe=list(probe=c("cg19403323"),geneNum=20),
             category="TN", dir.out ="./ELMER.example/Result/LUSC", save=FALSE)

## ----figure2, fig.height=4,fig.scap="Scatter Plot of One Pair", cache=TRUE----
scatter.plot(mee,byPair=list(probe=c("cg19403323"),gene=c("ID255928")),
             category="TN", save=FALSE,lm_line=TRUE)

## ----figure3, fig.height=4, fig.scap="TF expression vs. average DNA methylation", cache=TRUE----
load("ELMER.example/Result/LUSC/getMotif.hypo.enriched.motifs.rda")
scatter.plot(mee,byTF=list(TF=c("TP53","TP63","TP73"),
             probe=enriched.motif[["TP53"]]), category="TN",
             save=FALSE,lm_line=TRUE)

## ----makepair, cache=TRUE------------------------------------------------
# Make a "Pair" object for schematic.plot
pair <- fetch.pair(pair="./ELMER.example/Result/LUSC/getPair.hypo.pairs.significant.withmotif.csv",
                   probeInfo = "./ELMER.example/Result/LUSC/probeInfo_feature_distal.rda",
                   geneInfo = "./ELMER.example/Result/LUSC/geneAnnot.rda")

## ----figure4,fig.width=8, fig.height=3,fig.scap="Nearby Genes", cache=TRUE----
schematic.plot(pair=pair, byProbe="cg19403323",save=FALSE)

## ----figure5,fig.width=8, fig.height=3,fig.scap="Nearby Probes", cache=TRUE----
schematic.plot(pair=pair, byGene="ID255928",save=FALSE)

## ----figure6, fig.height=2,fig.scap="Motif Enrichment", cache=TRUE-------
motif.enrichment.plot(motif.enrichment="./ELMER.example/Result/LUSC/getMotif.hypo.motif.enrichment.csv",
                      significant=list(OR=1.3,lowerOR=1.3), dir.out ="ELMER.example/Result/LUSC",
                      label="hypo", save=FALSE)  ## different signficant cut off.

## ----figure7,fig.width=5.5, fig.height=3.5,fig.scap="TF Ranking", cache=TRUE----
load("./ELMER.example/Result/LUSC/getTF.hypo.TFs.with.motif.pvalue.rda")
TF.rank.plot(motif.pvalue=TF.meth.cor, motif="TP53", TF.label=list(TP53=c("TP53","TP63","TP73")),
            save=FALSE)

## ----finalsession, cache=TRUE--------------------------------------------
sessionInfo()

