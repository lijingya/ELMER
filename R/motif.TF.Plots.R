#'motif.enrichment.plot
#'@param motif.enrichment A data frame or file path of get.enriched.motif output 
#'motif.enrichment.csv file. See detial for the format of motif.enrichment if a data frame is specified.
#'@param significant A list to select subset of motif. Default is NULL. See detail
#'@param dir.out A path specify the directory to which the figures will be saved. 
#'Current directory is default.
#'@param save A logic. If true (default), figure will be saved to dir.out.
#'@param label A character labels the outputs figure.
#'@return A figure shows the enrichment level for selected motifs.
#'@details motif.enrichment If input data.frame object, it should contain "motif",
#' "OR", "lowerOR", "upperOR" columns. motif specifies name of motif; 
#' OR specifies Odds Ratio, lowerOR specifies lower boundary of OR (95%) ; 
#' upperOR specifies upper boundary of OR(95%).
#'@details significant A list used to select subset of motif.enrichment by the 
#'cutoff of OR, lowerOR, upperOR. significant=list(OR=1). More than one cutoff 
#'can be specified such as significant = list(OR=1, lowerOR=1,upperOR=4) 
#'@importFrom ggplot2 aes ggplot geom_point geom_errorbar coord_flip geom_abline
#'@export
#'@examples
#'motif.enrichment <- data.frame(motif=c("TP53","NR3C1","E2F1","EBF1","RFX5",
#'"ZNF143", "CTCF"),
#'OR=c(19.33,4.83,1, 4.18, 3.67,3.03,2.49),
#'lowerOR =c(10,3,1.09,1.9,1.5,1.5, 0.82),
#'upperOR =c(23,5,3,7,6,5,5),
#'stringsAsFactors=FALSE)
#'motif.enrichment.plot(motif.enrichment=motif.enrichment,
#'                      significant=list(OR=3),
#'                      label="hypo", save=FALSE)
motif.enrichment.plot <- function(motif.enrichment, significant=NULL, 
                                  dir.out ="./", save=TRUE,label=NULL){
  if(missing(motif.enrichment)) stop("motif.enrichment is missing.")
  if(is.character(motif.enrichment)){
    motif.enrichment <- read.csv(motif.enrichment, stringsAsFactors=FALSE)
  }
  if(!is.null(significant)){
    for(i in names(significant)){
      motif.enrichment <- motif.enrichment[motif.enrichment[,i] > significant[[i]],]
    }
  } 
  motif.enrichment <- motif.enrichment[order(motif.enrichment$OR,decreasing = TRUE),]
  motif.enrichment$motif <- factor(motif.enrichment$motif,
                                   levels=as.character(motif.enrichment$motif[nrow(motif.enrichment):1]))
  limits <- aes(ymax = upperOR, ymin=lowerOR)
  P <- ggplot(motif.enrichment, aes(x=motif, y=OR))+
    geom_point()+
    geom_errorbar(limits, width=0.3)+
    coord_flip()+
    geom_abline(intercept = 1, slope = 0, linetype = "3313")+
    theme_bw()+
    theme(panel.grid.major = element_blank())+
    labs(xlab="motifs", ylab="Odds Ratio")
  if(save) ggsave(filename = sprintf("%s/%s.motif.enrichment.pdf",dir.out,label),
                  useDingbats=FALSE, plot=P,width=4, 
                  height = 2*round(nrow(motif.enrichment)/10))
  return(P)
}


#'TF.rank.plot
#'@importFrom ggplot2 scale_color_manual geom_vline geom_text position_jitter 
#'@importFrom ggplot2 annotation_custom plot.margin unit ggplot_gtable ggplot_build
#'@importFrom grid textGrob linesGrob grid.draw
#'@importFrom gridExtra arrangeGrob
#'@param motif.pvalue A matrix or a path specifying location of "XXX.with.pvalue.rda" 
#'which is output of getTF. 
#'@param motif A vector of charactor specify the motif to plot
#'@param TF.label A list show the label for each motif. If TF.label is not specified, 
#'the motif relevant TF and top3 TF will be labeled.
#'@param dir.out A path specify the directory to which the figures will be saved. 
#'Current directory is default.
#'@param save A logic. If true (default), figure will be saved to dir.out.
#'@return A plot shows the score (-log(P value)) of association between TF
#'expression and DNA methylation of the certain motif sites.
#'@export
#'@examples
#'load(system.file("extdata","getTF.hypo.TFs.with.motif.pvalue.rda",package="ELMER"))
#'TF.rank.plot(motif.pvalue=TF.meth.cor, motif="TP53", TF.label=list(TP53=c("TP53","TP63","TP73")),
#'             save=FALSE)
TF.rank.plot <- function(motif.pvalue, motif, TF.label, dir.out="./", save=TRUE){
  if(missing(motif.pvalue)) stop("motif.pvalue should be specified.")
  if(missing(motif)) stop("Please specify which motif you want to plot.")
  if(is.character(motif.pvalue)){
    newenv <- new.env()
    load(motif.pvalue, envir=newenv)
    motif.pvalue <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  }
  if(missing(TF.label)){
    newenv <- new.env()
    data("motif.relavent.TFs",package = "ELMER.data",envir=newenv)
    motif.relavent.TFs <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
    TF.label <- motif.relavent.TFs[motif]
    specify <- "No"
  }else{
    specify <- "Yes"
  }
  significant <- floor(0.05*nrow(motif.pvalue))
  motif.pvalue <- -log10(motif.pvalue)
  
  Plots <- list()
  for(i in motif){
    df <- data.frame(pvalue=motif.pvalue[,i], Gene=rownames(motif.pvalue), stringAsFactors=FALSE)
    df <- df[order(df$pvalue, decreasing = TRUE),]
    df$rank <- 1:nrow(df)
    TF.sub <- TF.label[[i]]
    if(specify %in% "No"){
      TF.sub <- TF.sub[TF.sub %in% df$Gene[1:floor(0.05*nrow(df))]]
    }
    df$label <- "No"
    df$label[df$Gene %in% TF.sub|df$rank %in% 1:3] <- "Yes"
    df.label <- data.frame(pvalue = df$pvalue[df$label %in% "Yes"], 
                           text=as.character(df$Gene[df$label %in% "Yes"]), 
                           x=which(df$label %in% "Yes"), stringsAsFactors = FALSE)
    P <- ggplot(df, aes(x=rank, y=pvalue, color=factor(label, levels = c("Yes","No"))))+
      scale_color_manual(values = c("red","black"))+
      geom_vline(xintercept=significant, linetype = "3313")+
      geom_point()+
      # geom_text(aes(x=x+100, y=pvalue,label=text),df.label,color="black",
      #          position=position_jitter(height=0.5,width = 0.3), size=2)+
      #     geom_segment(aes(x = x, y = pvalue, xend = , yend = ),df.label,color="black")+
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      theme(legend.position="none")+
      labs(x="Rank", y="-log10 P value", title=i)
    
    
    TF.text <- list()
    delta <- ceiling(max(df$pvalue))/45
    for(i in df$Gene[df$label %in% "Yes"]){
      yposition <- df$pvalue[df$Gene %in% i]
      xposition <- df$rank[df$Gene %in% i]
      if(length(TF.text)>0){
        if(last.y-delta < yposition)
          yposition <- last.y-delta
      }
      # Create the text Grobs
      TF.text[[i]] = textGrob(i,gp=gpar(fontsize=10))
      # Draw the plot
      P = P + annotation_custom(grob = TF.text[[i]],  xmin = -370, xmax = -370, 
                                ymin = yposition, ymax = yposition) +
        annotation_custom(grob = linesGrob(), xmin = -270, xmax = xposition, 
                          ymin = yposition, ymax =df$pvalue[df$Gene %in% i] )
      last.y <- yposition
    }
    if(yposition < 0 ){
      y.margin <- abs(yposition)
    }else{
      y.margin <- 0
    }
    P <- P + theme(plot.margin = unit(c(y.margin, 0.5, 0, 2), "cm"))
    
    # Code to override clipping
    gt <- ggplot_gtable(ggplot_build(P))
    gt$layout$clip[gt$layout$name=="panel"] <- "off"
    
    if(save){
      pdf(sprintf("%s/%s.TFrankPlot.pdf",dir.out,i),
          useDingbats=FALSE, width=7, height = 6+y.margin)
      grid.draw(gt)
      invisible(dev.off())
    }else{
      Plots[[i]] <- arrangeGrob(gt)
    }
  }
  if(!save) return(Plots)
}




