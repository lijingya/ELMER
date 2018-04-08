#' scatter.plot to plot scatter plots between gene expression and DNA methylation.
#' @description 
#' scatter.plot is a function to plot various scatter plots between gene expression and 
#' DNA methylation. When byPair is specified, scatter plot for individual probe-gene pairs
#' will be generated. When byProbe is specified, scatter plots for one probes with nearby
#' 20 gene pairs will be generated. When byTF is specified, scatter plot for TF expression 
#' and average DNA methylation at certain motif sites will be generated.
#' @importFrom ggplot2 ggsave
#' @usage 
#' scatter.plot(data, 
#'              byPair = list(probe = c(), gene = c()),
#'              byProbe = list(probe = c(), numFlankingGenes = 20), 
#'              byTF = list(TF = c(), probe = c()),
#'              category = NULL,
#'              correlation = FALSE,
#'              width = 7,
#'              height = 6,
#'              dir.out = "./",
#'              save = TRUE, ...)
#' @param data A multiAssayExperiment with DNA methylation and Gene Expression data. 
#' See \code{\link{createMAE}} function.
#' @param byPair A list: byPair =list(probe=c(),gene=c()); probe contains a vector 
#'of probes' name and gene contains a vector of gene ID. The length of probe 
#'should be the same with length of gene. Output see numFlankingGenes
#'@param byProbe A list byProbe =list(probe=c(), geneNum=20); probe contains 
#'a vector of probes'name and geneNum specify the number of gene near the probes 
#'will ploted. 20 is default for numFlankingGenes Output see detail.
#'@param byTF A list byTF =list(TF=c(), probe=c()); TF contains a vector of TF's 
#'symbol and probe contains the a vector of probes' name. Output see detail.
#'@param category A vector labels subtype of samples or a character which is the 
#'column name in the colData(data) in the multiAssayExperiment object. Once specified, samples 
#'will label different color. The color can be customized by using color.value. 
#'@param dir.out A path specify the directory to which the figures will be saved. 
#'Current directory is default.
#'@param save A logic. If true, figure will be saved to dir.out.
#'@param height PDF height
#'@param width PDF width
#'@param correlation Add pearson correlation values to the plot
#'@param ... color.value, lm_line in scatter function
#'@details byPair The output will be scatter plot for individual pairs.
#'@details byProbe The output will be scatter plot for the probe and nearby genes.
#'@details byTF The output will be scatter plot for the TFs and the average 
#'DNA methylation at the probes set specified in byTF list.
#'@return Scatter plots.
#'@importFrom MultiAssayExperiment sampleMap
#'@export
#'@author Lijing Yao (maintainer: lijingya@usc.edu)
#'@examples
#' data <- ELMER:::getdata("elmer.data.example")
#' scatter.plot(data,
#'             byProbe=list(probe=c("cg19403323"),numFlankingGenes=20), 
#'             category="definition", save=FALSE)
#' scatter.plot(data,byProbe=list(probe=c("cg19403323"),numFlankingGenes=20), 
#'             category="definition", save=TRUE) ## save to pdf
#' # b. generate one probe-gene pair
#' scatter.plot(data,byPair=list(probe=c("cg19403323"),gene=c("ENSG00000143322")),
#'              category="definition", save=FALSE,lm_line=TRUE) 
scatter.plot <- function(data,
                         byPair = list(probe = c(),
                                       gene = c()),
                         byProbe = list(probe = c(),
                                        numFlankingGenes = 20),
                         byTF = list(TF = c(),
                                     probe = c()), 
                         category=NULL, 
                         correlation = FALSE,
                         width = 7,
                         height = 6,
                         dir.out = "./", 
                         save = TRUE, 
                         ...){
  simpleCap <- function(x) {
    s <- x
    paste(toupper(substring(s, first = 1, last = 1)), tolower(substring(s, 2)),
          sep = "", collapse = " ")
  }
    if(missing(data)) stop("A data object should be included.")
  
  if(!is.null(category) && length(category)==1) { 
    
    if(! category %in% colnames(colData(data))) 
      stop("Cateogry not found in the  phenotypic data (colData(data)) ")
    if(is.null(category)) stop("Please, set category argument")
    legend.title <- simpleCap(category)
    samples <- sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"]
    category <- colData(data)[samples,category]
    category <- sapply(category, simpleCap)
  }
  
  if(length(byPair$probe) != 0){
    
    if(length(byPair$probe) != length(byPair$gene))
      stop("In pairs, the length of probes should be the same with the length of genes.")

    pb <-  txtProgressBar(min = 0, max = length(byPair$gene), 
                          title = "creating images", 
                          style = 3, initial = 0, char = "=")
    
    for(i in 1:length(byPair$probe)){
      setTxtProgressBar(pb, i)
      probe <- byPair$probe[i]
      gene <- byPair$gene[i]
      symbol <- getSymbol(data,geneID=gene)
      P <- scatter(meth     = assay(getMet(data)[probe,]),
                   exp      = assay(getExp(data)[gene,] ),
                   category = category, 
                   legend.title = legend.title,
                   correlation = correlation,
                   xlab     = sprintf("DNA methyation at %s",probe), 
                   ylab     = sprintf("%s gene expression",symbol), 
                   title    = sprintf("%s_%s",probe,symbol),
                   ...)
      if(save) ggsave(filename = sprintf("%s/%s_%s_bypair.pdf",dir.out,probe,symbol),
                      plot = P,useDingbats=FALSE, width=width, height = height)
    }
    close(pb)  
    
  }
  if(length(byProbe$probe) != 0){
    nearGenes <- GetNearGenes(data    = data,
                              probes  = byProbe$probe,
                              numFlankingGenes = byProbe$numFlankingGenes)
    for(i in byProbe$probe){
      probe <- i
      gene <- nearGenes[[i]]$GeneID
      symbol <- getSymbol(data,geneID=gene)
      exp <- assay(getExp(data)[gene,])
      meth <- assay(getMet(data)[byProbe$probe,])
      rownames(exp) <- symbol
      P <- scatter(meth     = meth, 
                   exp      = exp,
                   category = category,
                   legend.title = legend.title,
                   xlab     = sprintf("DNA methyation at %s", probe), 
                   ylab     = sprintf("Gene expression"), 
                   title    = sprintf("%s nearby %s genes", probe, byProbe$numFlankingGenes),
                   ...)
      if(save) ggsave(filename = sprintf("%s/%s_byprobe.pdf", dir.out, probe),
                      plot = P, useDingbats=FALSE, width=width, height = height)
    }
  }
  
  if(length(byTF$TF)!=0){
    
    meth <- colMeans(assay(getMet(data)[byTF$probe,]),na.rm = TRUE)
    gene <- getGeneID(data,symbol = byTF$TF)

    # Our input might not be mapped, we need to verify it    
    found <- NULL
    if(any(is.na(gene))){
      found <- !is.na(gene)
      message("Gene not found: ", byTF$TF[!found])
      gene <- na.omit(gene) # rm the one not found
    }
    
    exp <- assay(getExp(data)[gene,])
    
    if(nrow(exp) > 1){
      if(!is.null(found)) {
        rownames(exp) <- byTF$TF[found]
      } else {
        rownames(exp) <- byTF$TF
      }
    }
    
    P <- scatter(meth     = meth, 
                 exp      = exp,
                 category = category,
                 legend.title = legend.title,
                 xlab     = "Avg DNA methyation", 
                 ylab     = sprintf("TF expression"), 
                 title    = "TF vs avg DNA methylation",
                 ...)
    
    if(save) ggsave(filename = sprintf("%s/%s_byTF.pdf",dir.out,paste(byTF$TF,collapse = "_")),
                    plot = P,useDingbats = FALSE, width = max(6,3*(length(byTF$TF)%%5)), 
                    height = max(4,3 * ceiling(length(byTF$TF)/5)))
  }
  return(P)
}


#'scatter
#'@importFrom reshape melt.data.frame
#'@import ggplot2
#'@param meth A vector of number.
#'@param exp A vector of number or matrix with sample in column and gene in rows.
#'@param category A vector of sample labels.
#'@param legend.title Plot legend title 
#'@param xlab A character specify the title of x axis.
#'@param ylab A character specify the title of y axis.
#'@param title A character specify the figure title.
#'@param correlation Show person correlation values 
#'@param color.value A vector specify the color of each category, such as 
#color.value=c("Experiment"="red","Control"="darkgreen")
#'@param lm_line A logic. If it is TRUE, regression line will be added to the graph.
#'@return A ggplot figure object
scatter <- function(meth, 
                    exp, 
                    legend.title = "Legend",
                    category=NULL, 
                    xlab=NULL, 
                    ylab=NULL,
                    title=NULL,
                    correlation = FALSE,
                    color.value=NULL,
                    lm_line=FALSE){
  
  if(is.null(category)) category <- rep(1,length(meth))
  
  if(!is.vector(exp)){
    exp <- as.data.frame(t(exp))
    GeneID <- colnames(exp)
    exp$meth <- as.vector(meth)
    exp$category <- category
    df <- melt.data.frame(exp, measure.vars = GeneID)
    P <- ggplot(df, aes(x = meth, y = value, color = factor(category))) +
      geom_point(size = 0.9) +
      facet_wrap(facets = ~ variable, ncol = 5) +
      scale_x_continuous(limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),  
            legend.position="bottom",
            legend.key = element_rect(colour = 'white'), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
      labs(x=xlab,y=ylab,title=title) + scale_colour_discrete(name=legend.title) + 
       guides(colour = guide_legend(override.aes = list(size=4),
                                    title.position="top", 
                                    title.hjust =0.5)) 
    if(!is.null(color.value)) P <- P + scale_colour_manual(values = color.value)
    if(lm_line) P <- P + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x,data=df)
    if(correlation){
      cor <- cor.test(x = as.numeric(meth), 
                      y = as.numeric(exp[,GeneID]),
                      method=c("pearson"))
      corval <- round(cor$estimate,digits = 2)
      pvalue <- round(cor$p.value,digits = 4)
      P <- P + annotate("text",x=0,y=0,hjust=.2,label=paste0("PCC: ",corval," / P-value: ",pvalue))
    }
  } else {
    df <- data.frame(meth=meth,exp=exp,category=category)
    if(length(unique(df$category))==1){
      P <- ggplot(df, aes(x = meth, y = exp))
    } else {
      P <- ggplot(df, aes(x = meth, y = exp, color = factor(category)))
    }
    P <- P + geom_point() +
      scale_x_continuous(limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position="bottom",
            legend.key = element_rect(colour = 'white'),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
      labs(x=xlab,y=ylab,title=title) +  
      scale_colour_discrete(name=legend.title) + 
      guides(colour = guide_legend(override.aes = list(size=4),
                                   title.position="top", 
                                   title.hjust =0.5))  +
      scale_fill_discrete(guide=FALSE)  + guides(fill=FALSE) 
    if(lm_line){
      #       P <- P+ geom_text(aes(x =0.8 , y = max(exp)-0.5, label = lm_eqn(df)),
      #parse = TRUE,colour = "black")+
      P <- P + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x,data=df)
    }

    
  }
  return(P)
}

