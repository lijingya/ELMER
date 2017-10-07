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
             font = list(size = 24), showlegend = showlegend,
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
      plotly::export(p, file = filename, vwidth = 992 * 2, vheight = 744 * 2)
      message("Saved as ", filename)
    } else {
      return(p)
    }
  })
}



#' Heatmp  of pairs gene and probes anti-correlated
#' @description 
#' Heatmp plot of pairs gene and probes anti-correlated
#' @param data A MultiAssayExperiment with a DNA methylation martrix or a DNA methylation matrix
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param pairs List of probe and pair genes
#' @param filename File names (.pdf) to save the file (i.e. "plot.pdf"). If NULL return plot.
#' @return A heatmap
#' @import ComplexHeatmap circlize
#' @export
#' @author Tiago Chedraoui Silva (tiagochst at gmail.com)
#' #' @examples
#' \dontrun{
#'   data <- ELMER:::getdata("elmer.data.example")
#'   group.col <- "subtype_Expression.Subtype"
#'   group1 <- "classical"
#'   group2 <- "secretory"
#'   pairs <- data.frame(Probe = c("cg15924102","cg19403323", "cg22396959"),
#'   GeneID = c("ENSG00000196878", "ENSG00000009790", "ENSG00000009790" ),
#'   Symbol = c("TRAF3IP3","LAMB3","LAMB3"),
#'   Raw.p =c(0.001,0.00001,0.001),
#'   Pe = c(0.001,0.00001,0.001))
#'   heatmapPairs(data, group.col,group1,group2,pairs,filename = "heatmap.pdf")
#'   }
heatmapPairs <- function(data, 
                         group.col, 
                         group1, 
                         group2, 
                         pairs, 
                         filename = NULL) {
  
  if(missing(data)) stop("Please set data argument")
  if(missing(group.col)) stop("Please set group.col argument")
  if(missing(group1)) stop("Please set group1 argument")
  if(missing(group2)) stop("Please set group2 argument")
  if(missing(pairs)) stop("Please set probe argument")
  
  data <- data[,colData(data)[,group.col] %in% c(group1, group2)]
  meth <- assay(getMet(data))[pairs$Probe,]
  exp <- assay(getExp(data))[pairs$GeneID,]
  order <- c(which(colData(data)[,group.col] == group1),which(colData(data)[,group.col] == group2))
  
  # Create color
  col <- rainbow(length(unique(colData(data)[,c(group.col)])))
  names(col) <- unique(colData(data)[,c(group.col)]) 
  col.list <- list()
  col.list[[group.col]] <- col
  
  # Annotation track
  ha = HeatmapAnnotation(df = colData(data)[,c(group.col),drop = F], 
                                         col = col.list)
  ha2 = HeatmapAnnotation(df = colData(data)[,c(group.col),drop = F], 
                                          show_legend = F,
                                          col = col.list)
  ht_list = 
   Heatmap(meth, name = "methylation", 
                            col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
                            column_names_gp = gpar(fontsize = 8),
                            show_column_names = F,
                            column_order = order, 
                            show_row_names = F,
                            cluster_columns = F,
                            top_annotation = ha,
                            column_title = "Methylation", column_title_gp = gpar(fontsize = 10), 
                            row_title_gp = gpar(fontsize = 10)) +
   Heatmap(t(apply(exp, 1, scale)), 
                            name = "expression (z-score)", 
                            col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                            top_annotation = ha2,
                            show_row_names = F,
                            column_order = order, 
                            cluster_columns = F,
                            column_names_gp = gpar(fontsize = 8), 
                            show_column_names = F,
                            column_title = "Expression", 
                            column_title_gp = gpar(fontsize = 10)) +
   ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                  heatmap_legend_labels_gp = gpar(fontsize = 8))
  if(is.null(filename)) return(draw(ht_list, newpage = TRUE, 
                                    column_title = "Correspondence between  probe DNA methylation and distal gene expression", 
                                    column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    annotation_legend_side = "bottom"))
  pdf(filename, width = 8, height = 8)
  
  draw(ht_list, newpage = TRUE, 
                       column_title = "Correspondence between probe DNA methylation and distal gene expression", 
                       column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
                       annotation_legend_side = "bottom")
  dev.off()
}