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
#' available columns using: colnames(MultiAssayExperiment::pData(data)).
#' @param group1 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param group2 A group from group.col. ELMER will run group1 vs group2. 
#' That means, if direction is hyper, get probes
#' hypermethylated in group 1 compared to group 2.
#' @param diff.dir A character can be "hypo" or "hyper", showing differential 
#'  methylation dirction.  It can be "hypo" which is only selecting hypomethylated probes; 
#' "hyper" which is only selecting hypermethylated probes; 
#' @param percentage A number ranges from 0 to 1 specifying the percentage of samples 
#' from group1 and group2 that are used to identify the differential methylation. 
#' Default is 0.2 because we did not expect all cases to be from a single molecular 
#' subtype.But, If you are working with molecular subtypes please set it to 1.
#' @param title plot title
#' @return Box plot
#' @importFrom plotly plot_ly layout
#' @importFrom dplyr top_n filter select %>%
#' @export
#' @author Tiago Chedraoui Silva (tiagochst at gmail.com)
#' @examples
#' data(elmer.data.example)
#' groupCol <- "subtype_Expression.Subtype"
#' group1 <- "classical"
#' group2 <- "secretory"
#' metBoxPlot(data,
#'            groupCol = groupCol,
#'            group1 = group1,
#'            group2 = group2, 
#'            probe ="cg17898069",
#'            percentage = 0.2, 
#'            diff.dir = "hypo")
metBoxPlot <- function(data, 
                       groupCol, 
                       group1, 
                       group2, 
                       probe, 
                       percentage = 0.2, 
                       diff.dir = "hypo",
                       title = "Title") {
  if(missing(data)) stop("Please set data argument")
  if(missing(groupCol)) stop("Please set data argument")
  if(missing(group1)) stop("Please set data argument")
  if(missing(group2)) stop("Please set data argument")
  if(missing(probe)) stop("Please set data argument")
  
  aux <- data.frame("group" = pData(data)[,groupCol],
                    "DNA methylation beta value" = assay(getMet(data))[probe,]) %>% 
    filter(group %in% c(group1,group2)) %>% droplevels
  
  if(diff.dir == "hyper") {
    val.used <- rbind(aux %>% filter(group %in% group2) %>% top_n(ceiling(sum(aux$group == group2) * .2)) %>% select(2),
                      aux %>% filter(group %in% group1) %>% top_n(ceiling(sum(aux$group == group1) * .2)) %>% select(2))
  } else {
    val.used <- rbind(aux %>% filter(group %in% group2) %>% top_n(-ceiling(sum(aux$group == group2) * .2)) %>% select(2),
                      aux %>% filter(group %in% group1) %>% top_n(-ceiling(sum(aux$group == group1) * .2)) %>% select(2))
  }
  aux$used <- ".all"
  aux.used <- aux[aux$DNA.methylation.beta.value %in% val.used$DNA.methylation.beta.value,]
  aux.used$used <- percentage
  
  aux <- rbind(aux,aux.used)
  aux$group <- paste0(aux$group, aux$used)
  suppressWarnings({
    p <- plot_ly(data = aux, 
                 x = ~group, 
                 y = ~DNA.methylation.beta.value,
                 color = ~used, 
                 type = "box", boxpoints = "all", jitter = 0.3,
                 pointpos = -1.8) %>% layout(title = title)
  })
  p
}