#' scatter.plot to plot scatter plots between gene expression and DNA methylation.
#' @description 
#' scatter.plot is a function to plot various scatter plots between gene expression and 
#' DNA methylation. When byPair is specified, scatter plot for individual probe-gene pairs
#' will be generated. When byProbe is specified, scatter plots for one probes with nearby
#' 20 gene pairs will be generated. When byTF is specified, scatter plot for TF expression 
#' and average DNA methylation at certain motif sites will be generated.
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. 
#' See \code{\link{createMAE}} function.
#' @param group.col A column defining the groups of the sample. You can view the 
#' available columns using: colnames(MultiAssayExperiment::colData(data)).
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param diff.dir A character can be "hypo" or "hyper", showing differential 
#'  methylation dirction.  It can be "hypo" which is only selecting hypomethylated probes; 
#' "hyper" which is only selecting hypermethylated probes; 
#' @param minSubgroupFrac A number ranges from 0 to 1 specifying the percentage of samples 
#' from group1 and group2 that are used to identify the differential methylation. 
#' Default is 0.2 because we did not expect all cases to be from a single molecular 
#' subtype.But, If you are working with molecular subtypes please set it to 1.
#' @param min.samples Minimun number of samples to use in the analysis. Default 5.
#' If you have 10 samples in one group, percentage is 0.2 this will give 2 samples 
#' in the lower quintile, but then 5 will be used.
#' @param title plot title
#' @param save Save plot as PNG
#' @param filename File names (.png) to save the file (i.e. "plot.png")
#' @param legend.col legend title
#' @param probe Character with probe name (i.e. "cg24517858")
#' @return Box plot
#' @importFrom plotly plot_ly layout
#' @importFrom dplyr top_n filter select %>%
#' @export
#' @author Tiago Chedraoui Silva (tiagochst at gmail.com)
#' @examples
#' \dontrun{
#'   data <- ELMER:::getdata("elmer.data.example")
#'   group.col <- "subtype_Expression.Subtype"
#'   group1 <- "classical"
#'   group2 <- "secretory"
#'   metBoxPlot(data,
#'              group.col = group.col,
#'              group1 = group1,
#'              group2 = group2, 
#'              probe ="cg17898069",
#'              minSubgroupFrac = 0.2, 
#'              diff.dir = "hypo")
#'}
metBoxPlot <- function(data, 
                       group.col, 
                       group1, 
                       group2, 
                       probe, 
                       min.samples = 5,
                       minSubgroupFrac = 0.2, 
                       diff.dir = "hypo",
                       legend.col = NULL,
                       title = NULL,
                       filename = NULL,
                       save = TRUE) {
  if(missing(data)) stop("Please set data argument")
  if(missing(group.col)) stop("Please set group.col argument")
  if(missing(group1)) stop("Please set group1 argument")
  if(missing(group2)) stop("Please set group2 argument")
  if(missing(probe)) stop("Please set probe argument")
  
  if(is.null(filename)){
    if(is.null(legend.col)) filename <- paste0(group.col,"_",probe)
    if(!is.null(legend.col)) filename <- paste0(group.col,"_",probe,"_",legend.col)
    filename <- paste0(gsub("\\.","_",filename),".png")
  }
  
  if(is.null(legend.col)) { 
    aux <- data.frame("group" = colData(data)[,group.col],
                      "DNA methylation beta value" = assay(getMet(data))[probe,]) %>% 
      filter(group %in% c(group1,group2)) %>% droplevels
    pos <- 2
    showlegend <- FALSE
  } else {
    aux <- data.frame("group" = colData(data)[,group.col], 
                      "legend" = colData(data)[,legend.col], 
                      "DNA methylation beta value" = assay(getMet(data))[probe,]) %>% 
      filter(group %in% c(group1,group2)) %>% droplevels
    if(any(is.na(aux$legend))) aux$legend[is.na(aux$legend)] <- "NA"
    pos <- 3
    showlegend <- TRUE
  }  
  if(diff.dir == "hyper") {
    val.used <- rbind(aux %>% filter(group %in% group2) %>% top_n(min(min.samples,ceiling(sum(aux$group == group2) * .2))) %>% select(pos),
                      aux %>% filter(group %in% group1) %>% top_n(min(min.samples,ceiling(sum(aux$group == group1) * .2))) %>% select(pos))
  } else {
    val.used <- rbind(aux %>% filter(group %in% group2) %>% top_n(-min(min.samples,ceiling(sum(aux$group == group2) * .2))) %>% select(pos),
                      aux %>% filter(group %in% group1) %>% top_n(-min(min.samples,ceiling(sum(aux$group == group1) * .2))) %>% select(pos))
  }
  aux$used <- " (100%)"
  aux.used <- aux[aux$DNA.methylation.beta.value %in% val.used$DNA.methylation.beta.value,]
  aux.used$used <- paste0(" (", minSubgroupFrac * 100, "%",")")
  aux <- rbind(aux,aux.used)
  aux$group <- paste0(aux$group, aux$used)
  if(is.null(title)) title <- paste0("Boxplot DNA methylation (probe ", probe,")")
  if(is.null(legend.col)) { 
    legend.col <- "used"
    legend.title <- ""
  } else {
    legend.title <- legend.col
    legend.col <- "legend"
  }
  
  suppressWarnings({
    p <- plot_ly(data = aux, 
                 x = ~group, 
                 y = ~DNA.methylation.beta.value,
                 color = ~eval(as.name(paste(legend.col))), 
                 type = "box", boxpoints = "all", jitter = 0.3,
                 pointpos = -1.8) %>% 
      plotly::add_annotations( text=legend.title, xref="paper", yref="paper",
                               x=1.02, xanchor="left",
                               y=0.8, yanchor="bottom",    # Same y as legend below
                               legendtitle=TRUE, showarrow=FALSE ) %>% 
      layout(title = title, boxmode = "group", xaxis = list(title = ""), 
             legend=list(y=0.8, yanchor="top" ), 
             font = list(size = 12), 
             showlegend = showlegend,
             margin = list(
               l = 100,
               r = 300,
               b = 100,
               t = 100,
               pad = 4
             ))
    if(save){
      if (!requireNamespace("webshot", quietly = TRUE)) {
        stop("webshot package is needed for this function to work. Please install it and run webshot::install_phantomjs()",
             call. = FALSE)
      }
      plotly::export(p, file = filename, vwidth = 992 , vheight = 744 )
      message("Saved as ", filename)
    } else {
      return(p)
    }
  })
}

#' Heatmap of pairs gene and probes anti-correlated
#' @description 
#' Heatmp plot of pairs gene and probes anti-correlated
#' @param data A MultiAssayExperiment with a DNA methylation SummarizedExperiment (all probes) and a gene Expression SummarizedExperiment.
#' @param group.col A column from the sample matrix from the MultiAssayExperiment object. Accessed with colData(mae)  
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param subset Subset MAE object to keep only groups compared ?
#' @param pairs List of probe and pair genes
#' @param annotation.col A vector of columns from the sample matrix from the MultiAssayExperiment object. Accessed with colData(mae) 
#' to be added as annotation to the heatmap.
#' @param met.metadata A vector of metdatada columns available in the DNA methylation GRanges to should be added to the heatmap. 
#' @param exp.metadata A vector of metdatada columns available in the Gene expression GRanges to should be added to the heatmap.
#' @param width Figure width
#' @param height Figure height
#' @param filename File names (.pdf) to save the file (i.e. "plot.pdf"). If NULL return plot.
#' @param cluster.within.groups Cluster columns based on the groups
#' @param plot.distNearestTSS Plot track with distNearestTSS ?
#' @return A heatmap
#' @import ComplexHeatmap circlize
#' @importFrom stats hclust dist
#' @importFrom grid unit.c grobWidth textGrob
#' @importFrom plyr ddply .
#' @importFrom GenomicRanges distanceToNearest
#' @importFrom grDevices png
#' @export
#' @author Tiago Chedraoui Silva (tiagochst at gmail.com)
#' @examples
#' \dontrun{
#'   data <- ELMER:::getdata("elmer.data.example")
#'   group.col <- "subtype_Expression.Subtype"
#'   group1 <- "classical"
#'   group2 <- "secretory"
#'   pairs <- data.frame(Probe = c("cg15924102","cg19403323", "cg22396959"),
#'                       GeneID = c("ENSG00000196878", "ENSG00000009790", "ENSG00000009790" ),
#'                       Symbol = c("TRAF3IP3","LAMB3","LAMB3"),
#'                       Distance = c(6017,168499,0),
#'                       Raw.p = c(0.001,0.00001,0.001),
#'                       Pe = c(0.001,0.00001,0.001))
#'  heatmapPairs(data = data, group.col = group.col,
#'               group1 = group1, group2 = group2,
#'               annotation.col = c("ethnicity","vital_status","age_at_diagnosis"),
#'               pairs, filename = "heatmap.pdf")
#'   }
heatmapPairs <- function(data, 
                         group.col, 
                         group1, 
                         group2, 
                         pairs, 
                         subset = FALSE,
                         cluster.within.groups = TRUE,
                         plot.distNearestTSS = FALSE,
                         annotation.col = NULL, 
                         met.metadata = NULL,
                         exp.metadata = NULL,
                         width = 10,
                         height = 7,
                         filename = NULL) {
  
  if(missing(data)) stop("Please set data argument")
  if(missing(group.col)) stop("Please set group.col argument")
  if(missing(pairs)) stop("Please set probe argument")
  
  if((!"distNearestTSS" %in% colnames(pairs)) & plot.distNearestTSS) {
    # For a given probe and gene find nearest TSS
    pairs <- addDistNearestTSS(data, pairs)
  }
  if(!missing(group1) & subset){
    message("Subsetting")
    data <- data[,colData(data)[,group.col] %in% c(group1, group2)]
  }
  
  # Remove pairs to be ploted if not found in the object
  pairs <- pairs[pairs$Probe %in% rownames(getMet(data)),]
  pairs <- pairs[pairs$GeneID %in% rownames(getExp(data)),]
  
  meth <- assay(getMet(data))[pairs$Probe,]
  exp <- assay(getExp(data))[pairs$GeneID,]
  
  order <- NULL
  cluster_columns <- TRUE
  if(cluster.within.groups){
    message("Ordering groups")
    cluster_columns <- FALSE
    order <- unlist(plyr::alply(unique(colData(data)[,group.col]), 1 , function(x) {
      if(is.na(x)){
        idx <- which(is.na(colData(data)[,group.col]))
      } else {
        idx <- which(colData(data)[,group.col] == x)
      }
      aux <- na.omit(meth[,idx])
      order <- t(aux) %>% dist %>% hclust(method = "average")
      as.numeric(idx[order$order])
    }
    ))
  }
  
  
  
  # Create color
  colors <- c("#6495ED", "#8B2323", "#458B74", "#7AC5CD",  
              "#4F4F4F", "#473C8B", "#00F5FF", "#CD6889", 
              "#B3EE3A", "#7B68EE", "#CDAF95", "#0F0F0F", "#FF7F00", 
              "#00008B", "#5F9EA0", "#F0FFFF", "#8B6969", "#9FB6CD", "#D02090",
              "#FFFF00", "#104E8B", "#B22222", "#B3EE3A", "#FF4500", "#4F94CD", 
              "#40E0D0", "#F5FFFA", "#8B3A62", "#FF3030", "#FFFFFF",
              "#191970","#BC8F8F","#778899","#2F4F4F",
              "#FFE4E1", "#F5F5DC")
  l <- length(unique(colData(data)[,c(group.col)]))
  l.all <-  l
  col <- colors[1:l]
  names(col) <- unique(colData(data)[,c(group.col)]) 
  names(col)[is.na(names(col))] <- "NA"
  col.list <- list()
  col.list[[group.col]] <- col
  
  for(i in annotation.col){
    l <- length(unique(colData(data)[,c(i)]))
    p <- l/length(na.omit(colData(data)[,c(i)]))
    if(p < 0.3 & !is.numeric(colData(data)[,c(i)])){
      if(l.all + l <= length(colors)) {
        col <- colors[(l.all+1):(l.all + l)]
        l.all <- l.all + l
      } else {
        col <- colors[c((l.all+1):length(colors),1 + (l.all + l)%%length(colors))]
        l.all <- (l.all + l)%%30
      }
      n <- unique(colData(data)[,c(i)]) 
      n[is.na(n)] <- "NA"
      names(col) <- n
      col.list[[i]] <- col
    } else {
      message("Considering variable ", i, " as numeric")
      suppressWarnings({
        nb <- as.numeric(colData(data)[,c(i)])
        colData(data)[,c(i)] <- nb
        if(!all(na.omit(nb) >=0)){
          col <- circlize::colorRamp2(c(min(nb,na.rm = T),
                                        (max(nb,na.rm = T) + min(nb,na.rm = T))/2,
                                        max(nb,na.rm = T)), c(colors[(l.all+1)],"white", colors[(l.all + 2)]))
          l.all <- l.all + 2
        } else {
          col.list[[i]] <- col
          col <- circlize::colorRamp2(c(min(nb,na.rm = T),max(nb,na.rm = T)), c("white", colors[(l.all+1):(l.all + 1)]))
          l.all <- l.all + 1
        }
        col.list[[i]] <- col
        
      })
    }
  }
  # Annotation track
  ha = HeatmapAnnotation(df = colData(data)[,c(group.col,annotation.col),drop = F], 
                         col = col.list,
                         show_annotation_name = TRUE,
                         annotation_name_side = "left",
                         annotation_name_gp = gpar(fontsize = 6))
  ha2 = HeatmapAnnotation(df = colData(data)[,c(group.col,annotation.col),drop = F], 
                          show_legend = F,
                          col = col.list)
  ht_list <- 
    Heatmap(meth, 
            name = "DNA methylation level", 
            col = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "gold")),
            column_names_gp = gpar(fontsize = 8),
            show_column_names = FALSE,
            column_order = order, 
            show_row_names = FALSE,
            use_raster = TRUE,
            raster_device = c("png"),
            raster_quality = 2,
            cluster_columns = cluster_columns,
            cluster_rows = TRUE,
            row_names_gp = gpar(fontsize = 6),
            top_annotation = ha,
            column_title = "DNA methylation", 
            column_title_gp = gpar(fontsize = 10), 
            row_title_gp = gpar(fontsize = 10)) 
  
  if(!is.null(met.metadata)) {
    for(i in met.metadata)
      ht_list <- ht_list + 
        Heatmap(values(getMet(data)[pairs$Probe,])[i],
                name = i, 
                use_raster = TRUE,
                raster_device = c("png"),
                raster_quality = 2,
                width = unit(5, "mm"),
                column_title = "", 
                show_column_names = F
        ) 
  }
  
  ht_list <-  ht_list +
    Heatmap(t(apply(exp, 1, scale)), 
            name = "Expression (z-score)", 
            col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
            top_annotation = ha2,
            show_row_names = FALSE,
            use_raster = TRUE,
            raster_device = c("png"),
            raster_quality = 2,
            column_order = order, 
            cluster_columns = cluster_columns,
            column_names_gp = gpar(fontsize = 8), 
            show_column_names = FALSE,
            column_title = "Expression", 
            column_title_gp = gpar(fontsize = 10)) 
  
  if(!is.null(exp.metadata)) {
    for(i in exp.metadata)
      ht_list <- ht_list + 
        Heatmap(values(getExp(data)[pairs$Probe,])[i],
                name = i, 
                use_raster = TRUE,
                raster_device = c("png"),
                raster_quality = 2,
                width = unit(5, "mm"),
                column_title = "", 
                show_column_names = FALSE
        ) 
  }
  if(plot.distNearestTSS){
    ht_list <-  ht_list +
      Heatmap(log10(pairs$distNearestTSS + 1),
              name = "log10(distNearestTSS + 1)", 
              width = unit(5, "mm"),
              column_title = "", 
              show_column_names = FALSE,
              use_raster = TRUE,
              raster_device = c("png"),
              raster_quality = 2,
              col = colorRamp2(c(0, 8), c("white", "orange")),
              heatmap_legend_param = list(at = log10(1 + c(0, 10, 100, 1000, 10000, 100000, 1000000,10000000,100000000)), 
                                          labels = c("0", "10bp", "100bp", "1kb", "10kb", "100kb", "1mb","10mb","100mb"))) 
  }
  ht_list <- ht_list +
    ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                  heatmap_legend_labels_gp = gpar(fontsize = 10))
  if(is.null(filename)) return(ht_list)
  padding = unit.c(unit(2, "mm"), grobWidth(textGrob(paste(rep("a",max(nchar(c(group.col,annotation.col)))/1.5), collapse = ""))) - unit(1, "cm"),
                   unit(c(2, 2), "mm"))
  if(grepl("\\.pdf",filename)) {
    message("Saving as PDF")
    pdf(filename, width = width, height = height)
  }
  if(grepl("\\.png",filename)) { 
    message("Saving as PNG")
    if(width < 100) width <- 1000
    if(height < 100) height <- 1000
    png(filename, width = width, height = height)
  }
  draw(ht_list, 
       #padding = padding, 
       newpage = TRUE, 
       column_title = "Correspondence between probe DNA methylation and distal gene expression", 
       column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
       annotation_legend_side = "right")
  dev.off()
}



#' @title Heatmap for correlation between probes DNA methylation and a single gene expression.
#' @description 
#' This heatmap will sort samples by their gene expression and show the DNA methylation levels of the paired probes to that gene.
#' If no pairs are given, nearest probes will be selected. 
#' To use this function you MAE object (input data) will need all probes and not only the distal ones.
#' This plot can be used to evaluate promoter, and intro, exons regions and closer distal probes of a gene to verify if their
#' DNA methylation level is affecting the gene expression
#' @param data A MultiAssayExperiment with a DNA methylation SummarizedExperiment (all probes) and a gene Expression SummarizedExperiment.
#' @param group.col A column from the sample matrix from the MultiAssayExperiment object. Accessed with colData(mae)  
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param pairs List of probe and pair genes
#' @param GeneSymbol Gene Symbol
#' @param annotation.col A vector of columns from the sample matrix from the MultiAssayExperiment object. Accessed with colData(mae) 
#' to be added as annotation to the heatmap
#' @param met.metadata A vector of metdatada columns available in the DNA methylation GRanges to should be added to the heatmap. 
#' @param exp.metadata A vector of metdatada columns available in the Gene expression GRanges to should be added to the heatmap.
#' @param scatter.plot Plot scatter plots
#' @param dir.out Where to save the plots
#' @param width Figure width
#' @param height Figure height
#' @param filename File names (.pdf) to save the file (i.e. "plot.pdf"). If NULL return plot.
#' @return A heatmap
#' @import ComplexHeatmap circlize
#' @importFrom stats hclust dist
#' @importFrom grid unit.c grobWidth textGrob
#' @importFrom plyr ddply .
#' @importFrom GenomicRanges distanceToNearest
#' @importFrom grDevices png
#' @export
#' @author Tiago Chedraoui Silva (tiagochst at gmail.com)
#' @examples
#' \dontrun{
#'   data <- ELMER:::getdata("elmer.data.example")
#'   group.col <- "subtype_Expression.Subtype"
#'   group1 <- "classical"
#'   group2 <- "secretory"
#'   pairs <- data.frame(Probe = c("cg15924102","cg19403323", "cg22396959"),
#'                       GeneID = c("ENSG00000196878", "ENSG00000009790", "ENSG00000009790" ),
#'                       Symbol = c("TRAF3IP3","LAMB3","LAMB3"),
#'                       Distance = c(6017,168499,0),
#'                       Raw.p = c(0.001,0.00001,0.001),
#'                       Pe = c(0.001,0.00001,0.001))
#'  heatmapGene(data = data, 
#'              group.col = group.col,
#'              group1 = group1, 
#'              group2 = group2,
#'              pairs = pairs, 
#'              GeneSymbol = "LAMB3",
#'              annotation.col = c("ethnicity","vital_status"),
#'              filename = "heatmap.pdf")
#'  \dontrun{            
#'      heatmapGene(data = data, 
#'                  group.col = group.col,
#'                  group1 = group1, 
#'                  group2 = group2,
#'                  GeneSymbol = "LAMB3",
#'                  annotation.col = c("ethnicity","vital_status"),
#'                  filename = "heatmap_closer_probes.pdf")
#'  }
#'}
heatmapGene <- function(data, 
                        group.col, 
                        group1, 
                        group2, 
                        pairs, 
                        GeneSymbol,
                        scatter.plot = FALSE,
                        annotation.col = NULL, 
                        met.metadata = NULL,
                        exp.metadata = NULL,
                        dir.out = ".",
                        width = 10,
                        height = 7,
                        filename = NULL) {
  
  if(missing(data)) stop("Please set data argument")
  if(missing(group.col)) stop("Please set group.col argument")
  if(missing(GeneSymbol)) stop("Please set GeneSymbol argument")
  dir.create(dir.out,showWarnings = FALSE, recursive = TRUE)
  # If pairs are missing we will select the probes that are closer to the gene
  if(missing(pairs)) {
    
    # Probe and gene info
    probes.info <- rowRanges(getMet(data)) # 450K and hg38
    gene.info <- rowRanges(getExp(data))
    
    gene.location <- gene.info[gene.info$external_gene_name == GeneSymbol,]
    # Get closer genes
    gene.follow <- gene.info[follow(gene.info[gene.info$external_gene_name == GeneSymbol,],gene.info,ignore.strand=T),]$external_gene_name
    gene.precede <- gene.info[precede(gene.info[gene.info$external_gene_name == GeneSymbol,],gene.info,ignore.strand=T),]$external_gene_name
    gene.gr <- gene.info[gene.info$external_gene_name %in% c(gene.follow,gene.precede,GeneSymbol),]
    # Get the regions of the 2 nearest genes, we will get all probes on those regions
    regions <- range(gene.gr)
    df <- data.frame(chrom= unique(as.character(seqnames(regions))), 
                     start= min(start(regions)), 
                     end = max(end(regions)))
    regions <- as(df, "GRanges")
    
    p <- names(sort(subsetByOverlaps(probes.info,regions,ignore.strand=TRUE)))
    if(length(p) == 0) stop("No probes close to the gene were found")
    neargenes <- ELMER::GetNearGenes(probes = p,
                                     data = data,
                                     numFlankingGenes = 10)
    pairs <- plyr::ldply(neargenes)
    pairs <- pairs[pairs$Symbol %in% GeneSymbol,]
    pairs <- addDistNearestTSS(data, pairs)
    p <- rev(names(sort(probes.info[p], ignore.strand=TRUE)))
    pairs <- pairs[match(p,pairs$Probe),]
    if(scatter.plot){
      for(p in pairs$Probe){
        scatter.plot(data = mae,
                     byPair = list(probe = p, 
                                   gene = gene.location$ensembl_gene_id), # TNFRSF11A
                     category = group.col, 
                     dir.out = dir.out, 
                     save = TRUE, 
                     correlation = TRUE,
                     lm_line = TRUE)
      }
    }
  }
  if(!GeneSymbol %in%  unique(pairs$Symbol)) stop("GeneID not in the pairs")
  pairs <- pairs[pairs$Symbol == GeneSymbol,]
  
  if(!"distNearestTSS" %in% colnames(pairs)) {
    # For a given probe and gene find nearest TSS
    pairs <- addDistNearestTSS(data, pairs)
  }
  if(!(missing(group1) & missing(group2))){
    data <- data[,colData(data)[,group.col] %in% c(group1, group2)]
  } 
  meth <- assay(getMet(data))[pairs$Probe,]
  exp <- assay(getExp(data))[unique(pairs$GeneID),,drop = FALSE]
  
  # Ordering the heatmap
  # Split data into the two groups and sort samples by the expression of the gene  
  if(!(missing(group1) & missing(group2))){
    idx1 <- which(colData(data)[,group.col] == group1)
    aux <- na.omit(exp[,idx1])
    dist1 <- sort(aux, index.return = T)$ix
    
    idx2 <- which(colData(data)[,group.col] == group2)
    aux <- na.omit(exp[,idx2])
    dist2 <- sort(aux, index.return = T)$ix
    
    order <- c(idx1[dist1],idx2[dist2])
  } else {
    aux <- na.omit(exp)
    dist1 <- sort(aux, index.return = T)$ix
    order <- dist1
  }  
  # Create color
  colors <- c("#6495ED", "#22b315", "#458B74", "#ffe300",  
              "#ff0000", "#473C8B", "#00F5FF", "#CD6889", 
              "#B3EE3A", "#7B68EE", "#CDAF95", "#0F0F0F", "#FF7F00", 
              "#00008B", "#5F9EA0", "#F0FFFF", "#8B6969", "#9FB6CD", "#D02090",
              "#FFFF00", "#104E8B", "#B22222", "#B3EE3A", "#FF4500", "#4F94CD", 
              "#40E0D0", "#F5FFFA", "#8B3A62", "#FF3030", "#FFFFFF")
  l <- length(unique(colData(data)[,c(group.col)]))
  l.all <-  l
  col <- colors[1:l]
  names(col) <- unique(colData(data)[,c(group.col)]) 
  col.list <- list()
  col.list[[group.col]] <- col
  
  for(i in annotation.col){
    l <- length(unique(colData(data)[,c(i)]))
    if(l < 10){
      if(l.all + l <= length(colors)) {
        col <- colors[(l.all+1):(l.all + l)]
        l.all <- l.all + l
      } else {
        col <- colors[c((l.all+1):length(colors),1 + (l.all + l)%%length(colors))]
        l.all <- (l.all + l)%%30
      }
      n <- unique(colData(data)[,c(i)]) 
      n[is.na(n)] <- "NA"
      names(col) <- n
      col.list[[i]] <- col
    } else {
      message("Considering variable ", i, " as numeric")
      suppressWarnings({
        nb <- as.numeric(colData(data)[,c(i)])
        colData(data)[,c(i)] <- nb
        if(!all(na.omit(nb) >=0)){
          col <- circlize::colorRamp2(c(min(nb,na.rm = T),
                                        (max(nb,na.rm = T) + min(nb,na.rm = T))/2,
                                        max(nb,na.rm = T)), c(colors[(l.all+1)],"white", colors[(l.all + 2)]))
          l.all <- l.all + 2
        } else {
          col.list[[i]] <- col
          col <- circlize::colorRamp2(c(min(nb,na.rm = T),max(nb,na.rm = T)), c("white", colors[(l.all+1):(l.all + 1)]))
          l.all <- l.all + 1
        }
        col.list[[i]] <- col
        
      })
    }
  }
  
  # Annotation track
  ha = HeatmapAnnotation(df = colData(data)[,c(group.col,annotation.col),drop = F], 
                         col = col.list,
                         show_annotation_name = TRUE,
                         annotation_name_side = "left",
                         annotation_height = unit(c(rep(0.5,length(annotation.col) + 1), 3), "cm"),
                         GeneExpression = anno_points(exp, size = unit(0.5, "mm"),axis = T, axis_side ="right"),
                         annotation_name_gp = gpar(fontsize = 6))
  
  bottom_annotation_height = unit(3, "cm")
  
  ha2 = HeatmapAnnotation(df = colData(data)[,c(group.col,annotation.col),drop = F], 
                          show_legend = F,
                          col = col.list)
  ht_list = 
    Heatmap(meth, 
            name = "DNA methylation level", 
            col = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "gold")),
            column_names_gp = gpar(fontsize = 5),
            show_column_names = F,
            column_order = order, 
            show_row_names = TRUE,
            cluster_columns = F,
            row_names_side = "left",
            cluster_rows = F,
            row_names_gp = gpar(fontsize = 5),
            top_annotation = ha,
            column_title_gp = gpar(fontsize = 10), 
            row_title_gp = gpar(fontsize = 10)) 
  
  if(!is.null(met.metadata)) {
    for(i in met.metadata)
      ht_list <- ht_list + 
        Heatmap(values(getMet(data)[pairs$Probe,])[i],
                name = i, 
                width = unit(5, "mm"),
                column_title = "", 
                show_column_names = F
        ) 
  }
  if(!is.null(exp.metadata)) {
    for(i in exp.metadata)
      ht_list <- ht_list + 
        Heatmap(values(getExp(data)[pairs$Probe,])[i],
                name = i, 
                width = unit(5, "mm"),
                column_title = "", 
                show_column_names = F
        ) 
  }
  
  ht_list <- ht_list +
    ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                  heatmap_legend_labels_gp = gpar(fontsize = 10))
  if(is.null(filename)) return(ht_list)
  padding = unit.c(unit(2, "mm"), grobWidth(textGrob(paste(rep("a",max(nchar(c(group.col,annotation.col)))/1.15), collapse = ""))) - unit(1, "cm"),
                   unit(c(2, 2), "mm"))
  if(grepl("\\.pdf",filename)) {
    message("Saving as PDF")
    pdf(sprintf("%s/%s",dir.out,filename), width = width, height = height)
  }
  if(grepl("\\.png",filename)) { 
    message("Saving as PNG")
    if(width < 100) width <- 1000
    if(height < 100) height <- 1000
    png(sprintf("%s/%s",dir.out,filename), width = width, height = height)
  }
  draw(ht_list, padding = padding, newpage = TRUE, 
       column_title = paste0("Correspondence between probe DNA methylation and ",  GeneSymbol," expression"), 
       column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
       annotation_legend_side = "right", 
       heatmap_legend_side = "bottom")
  dev.off()
}


#' @title Create a junction track for IGV visualization of interection
#' @description 
#' Create a junction track for IGV visualization of interection
#' @param pairs A data frame output from getPairs function
#' @param filename Filename (".bed")
#' @param met.platform DNA methyaltion platform to retrieve data from: EPIC or 450K (default)
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @param color.track A color for the track (i.e blue, red,#272E6A)
#' @param gene.symbol Filter pairs to a single gene.
#' @param track.name Track name
#' @param all.tss A logical. If TRUE it will link probes to all TSS  of a gene (transcript level), if FALSE 
#' it will link to the promoter region of a gene (gene level).
#' @importFrom dplyr %>%
#' @importFrom grDevices col2rgb
#' @importFrom utils write.table
#' @export
#' @author Tiago Chedraoui Silva (tiagochst at gmail.com)
#' @examples
#'  \dontrun{            
#' data <- ELMER:::getdata("elmer.data.example")
#' nearGenes <-GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
#'                          geneAnnot=getExp(data))
#' Hypo.pair <- get.pair(data=data,
#'                        nearGenes=nearGenes,
#'                        permu.size=5,
#'                        group.col = "definition",
#'                        group1 = "Primary solid Tumor", 
#'                        group2 = "Solid Tissue Normal",
#'                        raw.pvalue = 0.2,
#'                        Pe = 0.2,
#'                        dir.out="./",
#'                        label= "hypo")
#'  createIGVtrack(Hypo.pair,platform = "450K", genome = "hg38")
#'  }
createIGVtrack <- function(pairs,
                           met.platform = "450K",
                           genome = "hg38",
                           filename = "ELMER_interactions.bed",
                           color.track = "black",
                           track.name = "junctions",
                           gene.symbol = NULL,
                           all.tss = TRUE){
  if(all.tss){
    tss <- getTSS(genome = genome)
    tss <- tibble::as.tibble(tss)
  } else {
    tss <- TCGAbiolinks::get.GRCh.bioMart(genome = genome,as.granges = TRUE)
    tss <- tibble::as.tibble(promoters(tss,upstream = 0,downstream = 0))
    tss$transcription_start_site <- tss$start
  }
  
  if(!is.null(gene.symbol)) {
    if(!gene.symbol %in% pairs$Symbol) stop("Gene link with that gene symbol")
    pairs <- pairs[pairs$Symbol == gene.symbol,]
  }
  met.metadata <- getInfiniumAnnotation(plat = met.platform,genome = genome)
  met.metadata <- as.data.frame(met.metadata,row.names = names(met.metadata))
  met.metadata$Probe <- rownames(met.metadata)
  pairs <- merge(pairs, tss, by.x = "GeneID", by.y = "ensembl_gene_id",all.x = TRUE) %>% 
    merge(met.metadata,  by = "Probe", all.x = TRUE)
  
  pairs$ID <- paste0(pairs$Probe, ".", pairs$Symbol)
  pairs$geneCordinates <- paste0(0,",",pairs$width.x)
  pairs$strand <- "*"
  pairs$Raw.p <- as.integer(4)
  pairs$RGB <- paste(col2rgb(color.track)[,1],collapse = ",")
  pairs$block_counts <- 2
  pairs$block_sizes <- "0,0"
  # [seqname] [start] [end] [id] [score] [strand] 
  # [thickStart] [thickEnd] [r,g,b] [block_count] [block_sizes] [block_locations]
  
  pairs <- pairs[,c("seqnames.y", # [seqname]
                    "start.y",    # [start] # probe
                    "transcription_start_site",        # [end]
                    "ID",         # ID
                    "Raw.p",      # Depth
                    "strand",     # Strand
                    "start.x",    # thickStart
                    "end.x",      # thickEnd
                    "RGB",        # color 
                    "block_counts",    # block_count
                    "block_sizes",     # block_sizes
                    "block_counts")] # block_locations
  pairs <- na.omit(pairs)
  header <- paste0('track name="',track.name,'" graphType=junctions')
  unlink(filename,force = TRUE)
  cat(header,file = filename,sep="\n",append = TRUE)
  readr::write_delim(pairs, path = filename, append = TRUE)
}

#' @title Create a bigwig file for IGV visualization of DNA methylation data (Array)
#' @description 
#' Create a bigwig for IGV visualization of DNA methylation data (Array)
#' @param data A matrix
#' @param filename Filename (".bed")
#' @param met.platform DNA methyaltion platform to retrieve data from: EPIC or 450K (default)
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @importFrom plyr a_ply
#' @importFrom rtracklayer export.wig
#' @importFrom GenomeInfoDb Seqinfo
#' @export
#' @author Tiago Chedraoui Silva (tiagochst at gmail.com)
#' @examples
#'  \dontrun{            
#'  data <- assay(getMet(ELMER:::getdata("elmer.data.example")))
#'  createBigWigDNAmetArray(data = data, met.platform = "450K", genome = "hg38")
#'  }
createBigWigDNAmetArray <- function(data = NULL,
                                      genome = "hg38",
                                      met.platform = "450K",
                                      track.names = NULL,
                                      dir = "IGV_tracks"){
  
  # where we will save the several tracks
  dir.create(dir,recursive = TRUE, showWarnings = FALSE)
  
  # get genomic information for array
  message("Preparing array metadata")
  metadata <- getInfiniumAnnotation(plat = met.platform,genome = genome)
  metadata <- metadata[!metadata$MASK_general,]
  values(metadata) <- NULL
  metadata$score <- NA
  strand(metadata) <- "*"
  metadata <- keepStandardChromosomes(metadata, pruning.mode="coarse")
  seqinfo(metadata) <- Seqinfo(genome=genome) %>% keepStandardChromosomes
  
  message("Creating bigwig tracks")
  # for each samples create the track
  plyr::a_ply(colnames(data),1,.fun = function(sample){
    idx <- which(sample == colnames(data))
    met <- data[,sample]
    metadata <- metadata[names(metadata) %in% rownames(data),]
    data <- data[rownames(data) %in% names(metadata),]
    metadata$score <- met[match(rownames(data),names(metadata))]
    metadata <- metadata[!is.na(metadata$score),]
    if(!is.null(track.names)) {
      filename <- file.path(dir,track.names[idx])
    } else {
      filename <- file.path(dir,paste0(sample,".bw"))
    }
    message("\nSaving: ", filename)
    rtracklayer::export.bw(object = metadata,con = filename)
  },.progress = "time")
}




#' @title Creating matrix for MR TF heatmap    
#' @description Code used to create matrix for MR TF heatmap    
#' @param dir Vector ofr directory with results
#' @param classification Consider family or subfamily
#' @examples 
#' \dontrun{
#' elmer.results <- dirname(
#' dir(path = "analysis",
#'    pattern = "*.hypo.pairs.significant.csv",
#'    recursive = T,
#'    full.names = T,
#'    all.files = T))
#' tabs <- get.tabs(dir = elmer.results, classification = "subfamily")
#' }
get.tabs <- function(dir, classification = "family"){
  tab <- get.tab(dir,classification)
  tab.pval <- get.tab.pval(dir,classification,tab)
  tab.or <- get.tab.or(dir,classification,tab)
  tf.or.table <- get.tf.or.table(dir,classification,tab)
  return(list("tab" = tab,
              "tab.pval" = tab.pval,
              "tab.or" = tab.or,
              "tf.or.table" = tf.or.table))
}    

#' @importFrom dplyr full_join as_data_frame
#' @importFrom purrr reduce
get.tab <- function(dir,classification){
  message("o Creating TF binary matrix")
  tab <- purrr::reduce(lapply(dir, 
                              function(x) {
                                ret <- summarizeTF(path = x,
                                                   classification = classification,
                                                   top = TRUE)
                                colnames(ret)[2] <- x
                                return(ret)
                              }),
                       dplyr::full_join)
  tab[tab == "x"] <- 1
  tab[is.na(tab)] <- 0
  rownames(tab) <- tab$TF
  tab$TF <- NULL
  for(i in 1:ncol(tab)){
    tab[,i] <- as.numeric(tab[,i])
  }
  tab <- tab[rowSums(tab) > 0,]
  return(tab)
}    

#' @importFrom tidyr separate_rows gather
get.tab.or <- function(dir, classification, tab){
  col <- ifelse(classification == "family","top.potential.TF.family", "top.potential.TF.subfamily")
  message("o Creating TF OR matrix")
  # For each of those analysis get the enriched motifs and the MR TFs.
  tab.or <- plyr::adply(dir,
                        .margins = 1,
                        function(path){
                          TF <- readr::read_csv(dir(path = path, pattern = ".significant.TFs.with.motif.summary.csv",
                                                    recursive = T,full.names = T),col_types = readr::cols())
                          motif <- readr::read_csv(dir(path = path, pattern = ".motif.enrichment.csv",
                                                       recursive = T,full.names = T),col_types = readr::cols())
                          z <- tidyr::separate_rows(TF, col, convert = TRUE,sep = ";") %>% dplyr::full_join(motif)
                          z <- z[order(-z$OR),] # Drecreasing order of OR
                          OR <- z[match(rownames(tab),z[[col]]),] %>% pull(OR)
                          OR
                        },.id = NULL
  )
  tab.or <- t(tab.or)
  rownames(tab.or) <- rownames(tab)
  colnames(tab.or) <- colnames(tab)
  return(tab.or)
}

#' @importFrom tidyr separate_rows gather
get.tab.pval <- function(dir,classification,tab){ 
  col <- ifelse(classification == "family","top.potential.TF.family", "top.potential.TF.subfamily")
  message("o Creating TF FDR matrix")
  # For each of those analysis get the correlation pvalue of MR TFs exp  vs Avg DNA met. of inferred TFBS.
  tab.pval <- plyr::adply(dir,
                          .margins = 1,
                          function(path){
                            TF <- readr::read_csv(dir(path = path, 
                                                      pattern = ".significant.TFs.with.motif.summary.csv",
                                                      recursive = T,full.names = T),col_types = readr::cols())
                            motif <- readr::read_csv(dir(path = path, 
                                                         pattern = ".motif.enrichment.csv",
                                                         recursive = T,full.names = T),col_types = readr::cols())
                            # For each TF breaks into lines (P53;P63) will create two lines. Then merge with motif to get OR value
                            z <- tidyr::separate_rows(TF, col, convert = TRUE,sep = ";") %>% dplyr::full_join(motif)
                            z <- z[order(-z$OR),] # Drecreasing order of OR
                            colnames(z)[grep(paste0("^",col),colnames(z))] <- "TF"
                            motif <- z[match(rownames(tab),z$TF),] %>% pull(motif) # get higher OR for that TF binding
                            TF.meth.cor <- get(load(dir(path = path, pattern = ".TFs.with.motif.pvalue.rda", recursive = T, full.names = T)))
                            TF.meth.cor <- cbind(TF.meth.cor,TF = rownames(TF.meth.cor))
                            TF.meth.cor <- dplyr::as_data_frame(TF.meth.cor)
                            # Create table TF, motif, FDR
                            TF.meth.cor <- TF.meth.cor  %>% tidyr::gather(key = motif, value = FDR, -TF) 
                            # get FDR for TF and motif
                            FDR <- TF.meth.cor[match(paste0(rownames(tab),motif),paste0(TF.meth.cor$TF,TF.meth.cor$motif)),]  %>% dplyr::pull(FDR)
                            as.numeric(FDR)
                          },.id = NULL
  )
  
  tab.pval <- t(tab.pval)
  rownames(tab.pval) <- rownames(tab)
  colnames(tab.pval) <- colnames(tab)
  return(tab.pval)
}


#' @importFrom DelayedArray rowMins
get.top.tf.by.pval <- function(tab.pval,top = 5){
  labels <- c()
  for(i in 1:ncol(tab.pval)){
    labels <- sort(
      unique(
        c(labels,
          rownames(tab.pval)[head(sort(DelayedArray::rowMins(tab.pval[,i,drop = F],na.rm = T), 
                                       index.return = TRUE, 
                                       decreasing = F)$ix, 
                                  n = top)]
        )
      )
    )
  }
  
  return(labels)
}


get.tf.or.table <-function(dir,classification,tab){ 
  col <- ifelse(classification == "family","top.potential.TF.family", "top.potential.TF.subfamily")
  tf.or.table <- plyr::adply(dir,
                             .margins = 1,
                             function(path){
                               TF <- readr::read_csv(dir(path = path, 
                                                         pattern = ".significant.TFs.with.motif.summary.csv",
                                                         recursive = T,full.names = T),col_types = readr::cols())
                               motif <- readr::read_csv(dir(path = path, 
                                                            pattern = ".motif.enrichment.csv",
                                                            recursive = T,full.names = T),col_types = readr::cols())
                               
                               # For each TF breaks into lines (P53;P63) will create two lines. Then merge with motif to get OR value
                               z <- tidyr::separate_rows(TF, col, convert = TRUE,sep = ";") %>% dplyr::full_join(motif)
                               z <- z[order(-z$OR),] # Drecreasing order of OR
                               colnames(z)[grep(paste0("^",col),colnames(z))] <- "TF"
                               motif <- z[match(rownames(tab),z$TF),] %>% pull(motif) # get higher OR for that TF binding
                               TF.meth.cor <- get(load(dir(path = path, 
                                                           pattern = ".TFs.with.motif.pvalue.rda", 
                                                           recursive = T, full.names = T))
                               )
                               TF.meth.cor <- cbind(TF.meth.cor,TF =rownames(TF.meth.cor))
                               TF.meth.cor <- dplyr::as_data_frame(TF.meth.cor)
                               # Create table TF, motif, FDR
                               TF.meth.cor <- TF.meth.cor %>% tidyr::gather(key = motif, value = FDR, -TF)
                               TF.meth.cor$FDR <- as.numeric(TF.meth.cor$FDR)
                               TF.meth.cor <- TF.meth.cor[order(TF.meth.cor$FDR,decreasing = F),]
                               TF.meth.cor <- TF.meth.cor[
                                 match(paste0(rownames(tab),motif),
                                       paste0(TF.meth.cor$TF,TF.meth.cor$motif)),] %>% 
                                 na.omit
                               TF.meth.cor$analysis <- path
                               TF.meth.cor
                             },.id = NULL
  ) 
  tf.or.table <- tf.or.table[order(tf.or.table$FDR,decreasing = F),]
  return(tf.or.table)
}