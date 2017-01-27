#' schematic.plot to plot schematic plots showing the locations of genes and probes.
#' @description 
#' schematic.plot is a function to plot schematic plots showing the locations of genes and probes.
#' @usage 
#' schematic.plot(data,
#'                group.col = NULL,
#'                pair, 
#'                byProbe, 
#'                byGeneID, 
#'                byCoordinate=list(chr=c(), start=c(), end=c()),
#'                dir.out="./",
#'                save=TRUE,...)
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#'@param data A Multi Assay Experiment object with DNA methylation and 
#' gene expression Summarized Experiment objects
#' @param pair A data frame with two columns: Probe and Gene ID (Ensemble gene ID)
#' This is the ouput of get.pairfunction.
#' @param group.col A column defining the groups of the sample. You can view the 
#' available columns using: colnames(MultiAssayExperiment::pData(data)).
#' @param byProbe A vector of probe names.
#' @param byGeneID A vector of gene ID
#' @param byCoordinate A list contains chr, start and end. 
#'byCoordinate=list(chr=c(),start=c(),end=c()).
#' @param ... Parameters for GetNearGenes
#' @param dir.out A path specify the directory for outputs. Default is current directory
#' @param save A logic. If true, figures will be saved to dir.out.
#' @details
#'byProbes:
#'When a vector of probes' name are provided,
#'function will produce schematic plots for each individual probes.
#'The schematic plot contains probe, nearby 20 (or the number of gene user specified.)
#'genes and the significantly linked gene to the probe.
#'
#'byGene:
#'When a vector of gene ID are provided, function will produce schematic plots 
#'for each individual genes. The schematic plot contains the gene and all the 
#'significantly linked probes.
#'
#'byCoordinate:
#'When a genomic coordinate is provided, function will 
#'produce a schematic plot for this coordinate. The schematic plot contains 
#'all genes and significantly linked probes in the range and the significant links.
#' @export
#' @import Gviz lattice
#' @examples 
#' data(elmer.data.example)
#' pair <- data.frame(Probe = c("cg19403323","cg19403323", "cg26403223"),
#'                    GeneID = c("ENSG00000196878", "ENSG00000009790", "ENSG00000009790" ),
#'                    Symbol = c("TRAF3IP3","LAMB3","LAMB3"), 
#'                    Pe = c(0.001,0.00001,0.001))
#' schematic.plot(data, pair, byProbe = "cg19403323")
#' schematic.plot(data, pair, byGeneID = "ENSG00000009790")
#' schematic.plot(data, pair, byCoordinate = list(chr="chr1", start = 209000000, end = 209960000))
schematic.plot <- function(data,
                           group.col = NULL,
                           pair, 
                           byProbe, 
                           byGeneID, 
                           byCoordinate=list(chr=c(), start=c(), end=c()),
                           dir.out="./",
                           save=TRUE,...){
  # Begin of new schematic plot
  # For a probe get nearby genes
  if(!missing(byProbe)){
    args <- list(...)
    params <- args[names(args) %in% c("geneNum","cores")]
    nearGenes <- do.call(GetNearGenes,c(list(TRange=rowRanges(getMet(data))[byProbe,],
                                             geneAnnot=rowRanges(getExp(data))),params))
    for(probe in byProbe){
      significant <- pair[pair$Probe==probe,]
      gene.gr <- rowRanges(getExp(data))[nearGenes[[probe]]$GeneID,]
      probe.gr <- rowRanges(getMet(data))[unique(nearGenes[[probe]]$Target),]
      schematic(data = data, gene.gr, probe.gr, significant, label=sprintf("%s/%s.schematic.byProbe",dir.out,probe), save=save, group.col = group.col)
    }
  }
  if(!missing(byGeneID)){
    for(gene in byGeneID){
      
      significant <- pair[pair$GeneID==gene,]
      if(nrow(significant) == 0) {
        warning(paste0("Gene ", gene, " is not in pair list. We cannot plot it."))
        next
      }
      gene.gr <- rowRanges(getExp(data))[gene,]
      probe.gr <- rowRanges(getMet(data))[significant$Probe,]
      schematic(data = data,gene.gr, probe.gr, significant,label=sprintf("%s/%s.schematic.byGene",dir.out,gene), save=save, group.col = group.col)
    }
    
  }
  if(length(byCoordinate$chr)!=0){
    for(i in 1:length(byCoordinate$chr)){
      coordinate <- GRanges(seqnames = byCoordinate$chr[i], 
                            ranges = IRanges(byCoordinate$start[i],byCoordinate$end[i]))
      probe.gr  <- rowRanges(getMet(data))[pair$Probe,]
      probe.gr <-  probe.gr[queryHits(findOverlaps(probe.gr, coordinate)),]
      gene.gr <- rowRanges(getExp(data))[queryHits(findOverlaps(rowRanges(getExp(data)), coordinate)),]
      significant <- pair[pair$GeneID %in% names(gene.gr) & pair$Probe %in% names(probe.gr),]
      schematic(data = data,gene.gr, probe.gr, significant,
                    label=sprintf("%s/%s_%s_%s.schematic.byCoordinate",
                                  dir.out,byCoordinate$chr[i],byCoordinate$start[i],
                                  byCoordinate$end[i]), save=save, group.col = group.col)
    }
  }
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
                      group.col = NULL){
  
  if(save) pdf(paste0(label,".pdf"))
  chr <- as.character(seqnames(probe.gr))
  
  idxTrack <- IdeogramTrack(genome = metadata(data)$genome, chromosome = chr)
  axTrack <- GenomeAxisTrack()
  
  # We will find which is the significant pairs of genes
  fill <- rep("blue", length(values(gene.gr)$ensembl_gene_id))
 
  sig.colors <- rainbow(length(unique(significant$Probe)))
  for(i in seq_len(length(unique(significant$Probe)))) {
    fill[values(gene.gr)$ensembl_gene_id %in% significant[significant$Probe %in% unique(significant$Probe)[i],"GeneID"]] <- sig.colors[i]
  }
  genetrack <- GeneRegionTrack(gene.gr, name = "Gene", 
                               fill = fill, 
                               symbol = values(gene.gr)$external_gene_name,
                               shape = "arrow")
  
  details <- function(identifier, ...) {
    d <- data.frame(signal = assay(getMet(data))[identifier, ], group = pData(data)[,group.col])
    print(densityplot(~signal, group = group, data = d, auto.key = TRUE,
                      main = list(label = identifier, cex = 0.7),
                      scales = list(draw = FALSE, x = list(draw = TRUE)),
                      ylab = "", xlab = ""), newpage = FALSE,
          prefix = "plot")
  }
  if(!is.null(group.col)){
  deTrack <- AnnotationTrack(range = probe.gr, 
                             genome = metadata(data)$genome, 
                             showId = FALSE, 
                             groupAnnotation = "group",
                             chromosome = chr, 
                             fill = sig.colors, # We need to select the color based on the group
                             detailsConnector.col="red",
                             detailsBorder.col="red",
                             col.line="red",
                             detailsBorder.lwd=0,
                             id = names(probe.gr), 
                             name = "probe details",
                             stacking = "squish", 
                             fun = details)
  plotTracks(list(idxTrack,  axTrack, genetrack,deTrack),  
             background.title = "darkblue",
             detailsBorder.col = "white",
             sizes=c(1,1,1,8), 
             details.ratio = 1,
             details.size = 0.8,
             col = NULL,
             geneSymbols=TRUE)
  
  } else {
    plotTracks(list(idxTrack,  axTrack, genetrack),  
               background.title = "darkblue",
               detailsBorder.col = "white",
               sizes=c(1,1,1), 
               details.ratio = 1,
               details.size = 0.8,
               col = NULL,
               geneSymbols=TRUE)
    
  }
  if(save) dev.off()
}