#'schematicPlot
#'@import S4Vectors IRanges GenomicRanges
#'@param pair A pair object. All slots should be included
#'@param byProbe A vector of probe names.
#'@param byGene A vector of gene ID
#'@param byCoordinate A list contains chr, start and end. 
#'byCoordinate=list(chr=c(),start=c(),end=c()).
#'@param ... Parameters for GetNearGenes
#'@param dir.out A path specify the directory for outputs. Default is current directory
#'@param save A logic. If true, figure will be saved to dir.out.
#'@details byProbes When a vector of probes' name are provided, 
#'function will produce schematic plot for each individual probes. 
#'The schematic plot contains probe, nearby 20 (or the number of gene user specified.) 
#'genes and the significantly linked gene with the probe.
#'@details byGene When a vector of gene ID are provided, function will produce 
#'schematic plot for each individual genes. The schematic plot contains the gene 
#'and all the significantly linked probes.
#'@details byCoordinate When a genomic coordinate is provided, function will 
#'produce a schematic plot for this coordinate. The schematic plot contains 
#'all genes and significantly linked probes in the range and the significant links.
#'@return A schematic plot will be produced.
#'@export
#'@examples
#'library(grid)
#'load(system.file("extdata","mee.example.rda",package = "ELMER"))
#'nearGenes <-GetNearGenes(TRange=getProbeInfo(mee,probe=c("cg00329272","cg19403323")),
#'                         geneAnnot=getGeneInfo(mee))
#'Hypo.pair <-get.pair(mee=mee,probes=c("cg00329272","cg19403323"),
#'                     nearGenes=nearGenes,permu.size=5,Pe = 0.2,dir.out="./",
#'                     label= "hypo")
#'pair <- fetch.pair(pair=Hypo.pair,
#'                   probeInfo = getProbeInfo(mee),
#'                   geneInfo = getGeneInfo(mee))
#'# a. generate schematic plot of one probe with nearby 20 genes and label 
#'#the gene significantly linked with the probe.
#'grid.newpage()
#'schematic.plot(pair=pair, byProbe="cg19403323" ,save=FALSE)
#'#b. generate schematic plot of ont gene with the probe which the gene significanlty linked to.
#'grid.newpage()
#'schematic.plot(pair=pair, byGene="ID255928",save=FALSE)
schematic.plot <- function(pair, byProbe, byGene, 
                           byCoordinate=list(chr=c(), start=c(), end=c()),
                           dir.out="./",save=TRUE,...){
  if(missing(pair)) stop("Pair option is empty.")
  if(nrow(getPair(pair))==0 | length(getProbeInfo(pair))==0 | length(getGeneInfo(pair))==0) 
    stop("All slot should be included in pair object.")
  args <- list(...)
  if(!missing(byProbe)){
    params <- args[names(args) %in% c("geneNum","cores")]
    nearGenes <- do.call(GetNearGenes,c(list(TRange=getProbeInfo(pair,probe=byProbe),
                                             geneAnnot=getGeneInfo(pair)),params))
    for(i in byProbe){
      significant <- getPair(pair,probe=i)
      schematic(probe.range= getProbeInfo(pair,probe=i), 
                gene.range=getGeneInfo(pair, geneID=nearGenes[[i]]$GeneID),
                special=list(names= c(i,significant$Symbol), 
                             colors=c("blue",rep("red",length(significant$Symbol)))),
                label=sprintf("%s/%s.schematic.byProbe",dir.out,i), save=save)
    }
  }
  if(!missing(byGene)){
    for(i in byGene){
      significant <- getPair(pair,geneID=i)
      if(nrow(significant) ==0 ){
        warning(sprintf("gene %s don't have signficantly linked probes.",i))
      }else{
        schematic(probe.range= getProbeInfo(pair,probe=unique(significant$Probe)), 
                  gene.range=getGeneInfo(pair, geneID=i),
                  special=list(names= c(getSymbol(pair,geneID=i),significant$Probe), 
                               colors=c("red",rep("blue",length(significant$Probe)))),
                  label=sprintf("%s/%s.schematic.byGene",dir.out,i), save=save)
      }
    }
  }
  
  if(length(byCoordinate$chr)!=0){
    for(i in 1:length(byCoordinate$chr)){
      coordinate <- GRanges(seqnames = byCoordinate$chr[i], 
                            ranges =IRanges(byCoordinate$start[i],byCoordinate$end[i]))
      p.over <- findOverlaps(getProbeInfo(pair), coordinate)
      probe.range <- getProbeInfo(pair,probe=unique(getPair(pair)$Probe), range=coordinate) ## probe need to be significant
      gene.range <- getGeneInfo(pair, range=coordinate)   ## gene don't need
      significant <- getPair(pair,geneID=unique(gene.range$GENEID),probe=unique(gene.range$name))
      schematic(probe.range=probe.range , gene.range=gene.range,
                special=list(names= c(unique(significant$Symbol),unique(significant$Probe)), 
                             colors=c(rep("red",length(unique(significant$Symbol))),
                                      rep("blue",length(unique(significant$Probe))))),
                label=sprintf("%s/%s_%s_%s.schematic.byCoordinate",
                              dir.out,byCoordinate$chr[i],byCoordinate$start[i],
                              byCoordinate$end[i]), 
                interaction=list(probe=significant$Probe, gene= significant$Symbol),save=save)
    }
  }
}


#'schematicPlot
#'@param probe.range A GRanges object contains probes coordinates.
#'@param gene.range A GRanges object contains gene TSS coordinates.
#'@param special A list: special=list(names=c(), colors=c()) show the name of 
#'feature you want to highlight and specify the color respectively.
#'@param interaction A list: interaction=list(probe=c(),gene=c(),colors=c()) 
#'show the interacted features and specify the color respectively.
#'@param label A character labels the outputs figure.
#'@param save A logic. If true, figure will be saved to dir.out.
#'@return A schematic plot.
#'@details byProbes When a vector of probes' name are provided, function will 
#'produce schematic plot for each individual probes. The schematic plot contains 
#'probe, nearby 20 (or the number of gene user specified.) genes and the significantly
#' linked gene with the probe.
#'@details byGene When a vector of gene ID are provided, function will produce
#' schematic plot for each individual genes. The schematic plot contains the gene 
#' and all the significantly linked probes.
#'@details byCoordinate When a genomic coordinate is provided, function will 
#'produce a schematic plot for this coordinate. The schematic plot contains all 
#'genes and significantly linked probes in the range and the significant links.
#'@import IRanges GenomeInfoDb GenomicRanges
#'@importFrom grid grid.text grid.lines grid.curve arrow unit grid.circle viewport gpar
schematic <- function(probe.range, gene.range, special=list(names=c(),colors=c()),
                      interaction=list(probe=c(),gene=c(),colors=c()) ,label, save=TRUE){
  if(!unique(seqnames(probe.range)) %in% unique(seqnames(gene.range))) 
    stop("probe and gene should be in the same chromosome.")
  Sequences <- data.frame(name = gene.range$SYMBOL, 
                          position=start(gene.range),
                          Type=rep("arrow",length(gene.range)),
                          dir=as.vector(strand(gene.range)), 
                          stringsAsFactors = FALSE )
  Sequences <- rbind(Sequences, data.frame(name = probe.range$name, 
                                           position=start(probe.range), 
                                           Type=rep("circle",length(probe.range)), 
                                           dir=as.vector(strand(probe.range)),
                                           stringsAsFactors = FALSE ))
  Sequences <- do.call(rbind,lapply(split(Sequences,Sequences$name), 
                                    function(x){return(x[1,])}))
  Sequences <- Sequences[order(Sequences$position),]
  Sequences$colors <- rep("black", nrow(Sequences))
  Sequences$colors[match(special$names, Sequences$name)] <- special$colors
  Range <- range(Sequences$position)
  Sequences$dis <- 0.15 + 0.75*(Sequences$position - Sequences$position[1])/(Range[2]-Range[1])
  if(save) pdf(paste0(label,".pdf"))
  vp <- viewport(height=3, width=6)
  grid.text(paste0(unique(seqnames(gene.range)),":",Range[1],"-",Range[2]),0.5,0.9,gp=gpar(lwd=3))
  grid.lines(c(0.2,0.3),c(0.75,0.75),gp=gpar(lwd=3))
  grid.text(paste0(round((Range[2]-Range[1])/1000000,digits=2),"Mb"),0.25,0.78,gp=gpar(lwd=3))
  grid.lines(c(0.1,0.9),c(0.3,0.3),gp=gpar(lwd=3))
  
  ## plot gene
  sub.seq <- Sequences[Sequences$Type %in% "arrow" & Sequences$dir %in% "-",]
  if(nrow(sub.seq) >0){
    grid.curve(x1=sub.seq$dis, y1=0.3, x2=sub.seq$dis+(-0.025), y2=0.35,
               arrow=arrow(length = unit(0.01,"npc")),shape=0,
               gp=gpar(col=sub.seq$colors,lwd=3))
  }
  sub.seq <- Sequences[Sequences$Type %in% "arrow" & Sequences$dir %in% "+",]
  if(nrow(sub.seq) >0){
    grid.curve(x1=sub.seq$dis, y1=0.3, x2=sub.seq$dis+(0.025), y2=0.35, 
               arrow=arrow(length = unit(0.01,"npc")),shape=0,
               gp=gpar(col=sub.seq$colors,lwd=3),curvature = -1)
  }
  
  ## plot probe
  sub.seq <- Sequences[Sequences$Type %in% "circle",]
  grid.circle(x = sub.seq$dis, y = 0.3,r=0.01,
              gp=gpar(col=sub.seq$col, fill=sub.seq$col))
  
  ## label gene 
  sub.seq <- Sequences[Sequences$Type %in% "arrow" & 
                         !Sequences$colors %in% "black",]
  n <- 0.6
  for(i in 1:nrow(sub.seq)){
    x <- sub.seq[i,]
    grid.lines(x=c(x$dis,x$dis), y=c(0.35,n),gp=gpar(lwd=2))
    grid.text(x$name,x$dis,n+0.02,gp=gpar(lwd=3))
    n <- n-0.05
    if(n < 0.45) n<- 0.6
  }
  
  
  ## label probe
  sub.seq <- Sequences[Sequences$Type %in% "circle" &
                         !Sequences$colors %in% "black",]
  n <- 0.25
  for(i in 1:nrow(sub.seq)){
    x <- sub.seq[i,]
    if(i > 1){
      if(x$dis-sub.seq$dis[i-1] < 0.15 & n == 0.25){
        n <- 0.22
      }else if(x$dis-sub.seq$dis[i-1] < 0.15 & n != 0.22){
        n <- 0.25
      }
    }
    grid.text(x$name,x$dis,n+0.02,gp=gpar(lwd=3))
  }
  ## interaction
  if(length(interaction$probe)!=0 & length(interaction$gene)!=0){
    if(length(colors)==0) colors<- rep("black", length(Sequences))
    df <- data.frame(probe.seq=Sequences$dis[match(interaction$probe, Sequences$name)],
                     gene.seq=Sequences$dis[match(interaction$gene, Sequences$name)],
                     colors=interaction$colors)
    df$dir <- df$probe.seq-df$gene.seq
    
    ## curvature -1
    sub.df <- df[df$dir < 0,]
    grid.curve(x1=sub.df$probe.seq, y1=0.3, x2=sub.df$gene.seq, y2=0.35,
               gp=gpar(col=as.character(sub.df$colors),lwd=3),curvature = -1,
               ncp=8,square=0,angle=70)
    sub.df <- df[df$dir > 0,]
    grid.curve(x1=sub.df$probe.seq, y1=0.3, x2=sub.df$gene.seq, y2=0.35,
               gp=gpar(col=as.character(sub.df$colors),lwd=3),ncp=8,
               square=0,angle=-70)
  }
  if(save) dev.off()
  return(vp)
}
