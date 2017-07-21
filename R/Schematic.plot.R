#' schematic.plot to plot schematic plots showing the locations of genes and probes.
#' @description
#' schematic.plot is a function to plot schematic plots showing the locations of genes and probes.
#' @usage
#' schematic.plot(data,
#'                group.col = NULL,
#'                group1 = NULL,
#'                group2 = NULL,
#'                pair,
#'                byProbe,
#'                byGeneID,
#'                byCoordinate=list(chr=c(), start=c(), end=c()),
#'                statehub.tracks,
#'                dir.out="./",
#'                save=TRUE,...)
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#'@param data A Multi Assay Experiment object with DNA methylation and
#' gene expression Summarized Experiment objects
#' @param pair A data frame with three columns: Probe, Gene ID (Ensemble gene ID)
#' and Pe (empirical p-value). This is the ouput of get.pair function.
#' @param group.col A column defining the groups of the sample. You can view the
#' available columns using: colnames(MultiAssayExperiment::colData(data)).
#' @param group1 A group from group.col. ELMER will run group1 vs group2.
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2.
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.#' @param byProbe A vector of probe names.
#' @param byGeneID A vector of gene ID
#' @param byProbe A vector of probe names
#' @param byCoordinate A list contains chr, start and end.
#'byCoordinate=list(chr=c(),start=c(),end=c()).
#' @param ... Parameters for GetNearGenes
#' @param dir.out A path specify the directory for outputs. Default is current directory
#' @param statehub.tracks Relative path to a statehub track.
#' @param save A logic. If true, figures will be saved to dir.out.
#' @details
#' byProbes:
#'  When a vector of probes' name are provided,
#'  function will produce schematic plots for each individual probes.
#'  The schematic plot contains probe, nearby 20 (or the number of gene user specified.)
#'  genes and the significantly linked gene to the probe.
#'
#' byGene:
#'  When a vector of gene ID are provided, function will produce schematic plots
#'  for each individual genes. The schematic plot contains the gene and all the
#'  significantly linked probes.
#'
#' byCoordinate:
#'  When a genomic coordinate is provided, function will
#'  produce a schematic plot for this coordinate. The schematic plot contains
#'  all genes and significantly linked probes in the range and the significant links.
#' @export
#' @import Gviz lattice
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom  methods as is
#' @examples
#' data <- ELMER:::getdata("elmer.data.example")
#' pair <- data.frame(Probe = c("cg19403323","cg19403323", "cg26403223"),
#'                    GeneID = c("ENSG00000196878", "ENSG00000009790", "ENSG00000009790" ),
#'                    Symbol = c("TRAF3IP3","LAMB3","LAMB3"),
#'                    Raw.p =c(0.001,0.00001,0.001),
#'                    Pe = c(0.001,0.00001,0.001))
#' schematic.plot(data,
#'                group.col = "definition",
#'                group1 = "Primary solid Tumor",
#'                group2 = "Solid Tissue Normal",
#'                pair = pair,
#'                byProbe = "cg19403323")
#' schematic.plot(data,
#'                group.col = "definition",
#'                group1 = "Primary solid Tumor",
#'                group2 = "Solid Tissue Normal",
#'                pair = pair,
#'                byGeneID = "ENSG00000009790")
#'
#' schematic.plot(data,
#'                group.col = "definition",
#'                group1 = "Primary solid Tumor",
#'                group2 = "Solid Tissue Normal",
#'                pair = pair,
#'                byCoordinate = list(chr="chr1", start = 209000000, end = 209960000))
#' \dontrun{
#'    schematic.plot(data,
#'                   group.col = "definition",
#'                   group1 = "Primary solid Tumor",
#'                   group2 = "Solid Tissue Normal",
#'                   pair = pair,
#'                   byProbe = "cg19403323",
#'                   statehub.tracks = "hg38/ENCODE/mcf-7.16mark.segmentation.bed")
#' }
schematic.plot <- function(data,
                           group.col = NULL,
                           group1 = NULL,
                           group2 = NULL,
                           pair,
                           byProbe,
                           byGeneID,
                           byCoordinate=list(chr=c(), start=c(), end=c()),
                           statehub.tracks = NULL,
                           dir.out="./",
                           save=TRUE,...){
  # Begin of new schematic plot
  # For a probe get nearby genes
  args <- list(...)
  params <- args[names(args) %in% c("numFlankingGenes","cores")]
  
  extra.tracks <- args[grepl("track",names(args))]
  
  suppressWarnings({
    if(!missing(byProbe)){
      nearGenes <- do.call(GetNearGenes,c(list(TRange=rowRanges(getMet(data))[byProbe,],
                                               geneAnnot=rowRanges(getExp(data))),params))
      pb <-  txtProgressBar(min = 0, max = length(byProbe), title = "creating images", 
                            style = 3, initial = 0, char = "=")
      progress <- 0
      for(probe in byProbe){
        progress <- progress + 1
        setTxtProgressBar(pb, progress)
        significant <- pair[pair$Probe == probe,]
        gene.gr <- rowRanges(getExp(data))[nearGenes[[probe]]$GeneID,]
        probe.gr <- rowRanges(getMet(data))[unique(nearGenes[[probe]]$Target),]
        schematic(data = data,
                  gene.gr   = gene.gr,
                  probe.gr  =  probe.gr,
                  significant = significant,
                  label     = sprintf("%s/%s.schematic.byProbe",dir.out,probe),
                  statehub.tracks = statehub.tracks,
                  save      = save,
                  group.col = group.col,
                  group1    = group1,
                  group2    = group2,
                  extra.tracks = extra.tracks)
      }
      close(pb)  
    }
    if(!missing(byGeneID)){
      pb <-  txtProgressBar(min = 0, max = length(byGeneID), title = "creating images", 
                            style = 3, initial = 0, char = "=")
      progress <- 0
      for(gene in byGeneID){
        progress <- progress + 1
        setTxtProgressBar(pb, progress)
        significant <- pair[pair$GeneID==gene,]
        if(nrow(significant) == 0) {
          warning(paste0("Gene ", gene, " is not in pair list. We cannot plot it."))
          next
        }
        gene.gr <- rowRanges(getExp(data))[gene,]
        probe.gr <- rowRanges(getMet(data))[significant$Probe,]
        schematic(data = data,
                  gene.gr,
                  probe.gr,
                  significant,
                  label = sprintf("%s/%s.schematic.byGene",dir.out,gene),
                  save  = save,
                  statehub.tracks = statehub.tracks,
                  group.col = group.col,
                  group1 = group1,
                  group2 = group2,
                  extra.tracks = extra.tracks)
      }
      close(pb)
    }
    if(length(byCoordinate$chr)!=0){
      for(i in 1:length(byCoordinate$chr)){
        coordinate <- GRanges(seqnames = byCoordinate$chr[i],
                              ranges = IRanges(byCoordinate$start[i],byCoordinate$end[i]))
        probe.gr  <- rowRanges(getMet(data))[unique(pair$Probe),]
        probe.gr <-  probe.gr[queryHits(findOverlaps(probe.gr, coordinate)),]
        if(length(probe.gr) == 0) stop("No probes in that region")
        gene.gr <- rowRanges(getExp(data))[queryHits(findOverlaps(rowRanges(getExp(data)), coordinate)),]
        if(length(gene.gr) == 0) stop("No genes in that region")
        significant <- pair[pair$GeneID %in% names(gene.gr) & pair$Probe %in% names(probe.gr),]
        schematic(data = data,
                  gene.gr,
                  probe.gr,
                  significant,
                  label=sprintf("%s/%s_%s_%s.schematic.byCoordinate",
                                dir.out,byCoordinate$chr[i],byCoordinate$start[i],
                                byCoordinate$end[i]), save=save,
                  statehub.tracks = statehub.tracks,
                  group.col = group.col,
                  group1 = group1,
                  group2 = group2,
                  extra.tracks = extra.tracks)
      }
    }
  })
}

#' @importFrom grDevices rainbow
#' @importFrom GenomicRanges seqnames
#' @importFrom MultiAssayExperiment metadata
schematic <- function(data,
                      gene.gr,
                      probe.gr,
                      significant,
                      label,
                      save=TRUE,
                      statehub.tracks = NULL,
                      group.col = NULL,
                      group1 = NULL,
                      group2 = NULL,
                      extra.tracks = NULL){
  options(ucscChromosomeNames=FALSE)
  
  chr <- as.character(seqnames(probe.gr))
  
  idxTrack <- IdeogramTrack(genome = metadata(data)$genome, chromosome = chr)
  axTrack <- GenomeAxisTrack()
  
  # We will find which is the significant pairs of genes
  fill <- rep("blue", length(values(gene.gr)$ensembl_gene_id))
  
  for(i in seq_len(length(unique(significant$Probe)))) {
    fill[values(gene.gr)$ensembl_gene_id %in% significant[significant$Probe %in% unique(significant$Probe)[i],]$GeneID] <- "red"
  }
  genetrack <- GeneRegionTrack(gene.gr, 
                               name = "Gene",
                               fill = fill,
                               symbol = values(gene.gr)$external_gene_name,
                               shape = "arrow")
  
  details <- function(identifier, ...) {
    d <- data.frame(signal = assay(getMet(data))[identifier, ], group = colData(data)[,group.col])
    print( 
      bwplot(signal~group,
             data=d,
             xlab=group.col, ylab='DNA methylation levels',
             horizontal=FALSE,
             panel = function(..., box.ratio) {
               panel.violin(..., col = "lightblue",
                            varwidth = FALSE, box.ratio = box.ratio)
               panel.bwplot(..., col='black',
                            cex=0.8, pch='|', fill='gray', box.ratio = .1)
             },
             par.settings = list(box.rectangle=list(col='black'),
                                 plot.symbol = list(pch='.', cex = 0.1)),
             scales=list(x=list(rot=0, cex=0.5))),
      #densityplot(~signal, group = group, data = d, auto.key = TRUE,
      #                main = list(label = identifier, cex = 0.7),
      #                scales = list(draw = FALSE, x = list(draw = TRUE)),
      #                ylab = "", xlab = ""), 
      newpage = FALSE,
      prefix = "plot")
  }
  interactions.track <- c()
  if(nrow(significant) > 0 ) {
    if (!requireNamespace("GenomicInteractions", quietly = TRUE)) {
      stop("GenomicInteractions package is needed for this function to work. Please install it.",
           call. = FALSE)
    }
    genes.plot <- gene.gr[match(significant$GeneID,names(gene.gr))]
    genes.plot <- SummarizedExperiment::resize(genes.plot,width = 1)
    interactions <- GenomicInteractions::GenomicInteractions(genes.plot,
                                        probe.gr[match(significant$Probe,names(probe.gr))],
                                        experiment_name="Putative pair genes ",
                                        description="this is a test", counts=-log10(significant$Raw.p))
    interactions.track <-  GenomicInteractions::InteractionTrack(name="Putative pair genes\n (-log10 raw p-value)", interactions, chromosome=chr)
    displayPars(interactions.track) = list(col.interactions="red", 
                                           col.anchors.fill ="transparent",
                                           col.anchors.line = "transparent",
                                           interaction.dimension="height", 
                                           interaction.measure ="counts",
                                           anchor.height = 0.1)
  }
  probe.col <- "black"
  # StateHub tracks
  state.tracks <- c()
  if(!is.null(statehub.tracks)){
    base <- "http://s3-us-west-2.amazonaws.com/statehub-trackhub/tracks/5813b67f46e0fb06b493ceb0/"
    for(state in statehub.tracks){
      message("Adding stateHub track: ", state)
      bed <- paste0(base,state)
      if(!file.exists(basename(bed))) downloader::download(bed,basename(bed))
      if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("rtracklayer package is needed for this function to work. Please install it.",
             call. = FALSE)
      }
      state.chr <- rtracklayer::import.bed(basename(bed))
      state.chr <- state.chr[seqnames(state.chr) == chr]
      #state.chr <-  state[seqnames(state) == chr &
      #                      start(state) >= min(start(gene.gr) , start(probe.gr) ) &
      #                      end(state) <= max(end(gene.gr) , end(probe.gr) )]
      
      tracks <- plyr::alply(unique(state.chr$name), 1, function(x){
        aux <- state.chr[state.chr$name == x]
        AnnotationTrack(aux,name = paste0(state.chr@trackLine@name, "\n",x), 
                        stacking = "dense",
                        col = NULL,
                        col.line = NULL,
                        shape = "box",
                        fill = unique(aux$itemRgb))
      })
    }
    #all.states <- AnnotationTrack(state.chr,fill = state.chr$itemRgb, stacking = "squish")
    state.tracks <- c(state.tracks,tracks)
    message("Probes overlapping")
    print(IRanges::subsetByOverlaps(state.chr,probe.gr)$name)
  }
  if(save) pdf(paste0(label,".pdf"), height = max(5, 5 + rep(2,!is.null(group.col)),
                                                  floor(length(state.tracks)/2 + 5 + rep(2,!is.null(group.col)))))
  
  
  if(!is.null(group.col)){
    
    if(!is.null(group1) & !is.null(group1))
      data <- data[,colData(data)[,group.col] %in% c( group1, group2)]
    
    deTrack <- AnnotationTrack(range = probe.gr,
                               genome = metadata(data)$genome,
                               showId = FALSE,
                               groupAnnotation = "group",
                               chromosome = chr,
                               fill = probe.col,
                               detailsConnector.col="grey",
                               detailsBorder.col="grey",
                               col.line=probe.col,
                               detailsBorder.lwd=0,
                               detailsConnector.pch = NA,
                               id = names(probe.gr),
                               name = "Probe details",
                               stacking = "squish",
                               fun = details)
    plotTracks(c(list(idxTrack,  axTrack, deTrack), interactions.track,list(genetrack),extra.tracks,state.tracks),
               background.title = "darkblue",
               detailsBorder.col = "white",
               from = min(start(gene.gr) , start(probe.gr), end(gene.gr) , end(probe.gr)),
               to = max(start(gene.gr) , start(probe.gr), end(gene.gr) , end(probe.gr)),
               sizes=c(1,1,8,rep(2,length(interactions.track)),3,rep(2,length(extra.tracks)),rep(0.5,length(state.tracks))),
               extend.right = 10000,
               extend.left = 100000,
               details.ratio = 1,
               details.size = 0.9,
               baseline=0, innerMargin=0,
               col = NULL,
               #fontsize = 8,
               showBandId = TRUE, cex.bands = 0.5,
               title.width = 2,
               cex.title = 0.5,
               rotation.title=360,
               geneSymbols=TRUE)
  } else {
    atrack <- AnnotationTrack(probe.gr, name = "Probes",
                              genome = metadata(data)$genome,
                              chromosome = chr)
    plotTracks(c(list(idxTrack,  axTrack), interactions.track,list(atrack , genetrack),extra.tracks,state.tracks),
               background.title = "darkblue",
               from = min(start(gene.gr) , start(probe.gr), end(gene.gr) , end(probe.gr)),
               to = max(start(gene.gr) , start(probe.gr), end(gene.gr) , end(probe.gr)),
               extend.right = 10000,
               extend.left = 100000,
               baseline=0, innerMargin=0,
               showBandId = TRUE, cex.bands = 0.5,
               detailsBorder.col = "white",
               sizes=c(0.5,0.5,rep(2,length(interactions.track)),1,3,rep(2,length(extra.tracks)),rep(0.5,length(state.tracks))),
               details.ratio = 1,
               #fontsize = 8,
               rotation.title=360,
               title.width = 2,
               details.size = 0.8,
               col = NULL,
               cex.title = 0.5,
               geneSymbols=TRUE)
  }
  if(save) {
    message("Saving as: ", label)
    dev.off()
  }
}
