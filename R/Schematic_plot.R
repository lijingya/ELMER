###schematic figure funciton
#scale # per 1 unit
library(grid)
#' Generate random loci of genome.
#' @param target.range A GRange object showing target coordinate.
#' @param Gene.range A GRange oject contains coordinate of a list of genes to show on figure.
#' @param special A list. Symbols is specific the gene that you want to emphasize and colors to specific the color to each gene you want to emphasize.
#' @param target.col A color for target
#' @param save Specific save the graph object or not
#' @return A graph object or save as a pdf.

schematic <- function(target.range,Gene.range,special=list(Symbols=c(),colors=c()),target.col,fn,save=T){
  tmp <- Gene.range
  strand(tmp) <- "*"
  tmp <- sort(tmp)
  position <- follow(target.range,tmp)
  n <- 1
  Righttmp <- tmp[(position+1)]
  while(position+n+1 < length(tmp)+1){
    if(!as.character(tmp$SYMBOL[position+n+1]) %in% as.character(Righttmp$SYMBOL)){
      Righttmp <- c(Righttmp,tmp[position+n+1])
    }
    n <- n+1
  }
  if(length(as.character(Righttmp$SYMBOL)) < 20){
    n <-0
    while(as.character(tmp$SYMBOL[position-n]) %in% as.character(Righttmp$SYMBOL)){
      n <- n+1
    }
    Lefttmp <- tmp[position-n]
    while(position-n >= 1){
      if(!as.character(tmp$SYMBOL[position-n]) %in% as.character(Lefttmp$SYMBOL)){
        Lefttmp <- c(tmp[position-n],Lefttmp)
      }
      n <- n+1
    }
    Unique <- c(Lefttmp,Righttmp)
  }else{
    Unique <- Righttmp
  }
  
  if(length(Righttmp)==20){
    Range <- distance(target,Unique[length(Unique)])
    targetD <- 0.25 
  }else{
    Range <- distance(Unique[1],Unique[length(Unique)])
    targetD <- 0.25 + 0.7* distance(Unique[1],target.range)/Range
  }
  Distance <- 0.25 + 0.7* unlist(apply(matrix(1:length(Unique),ncol=1),1, function(x){ distance(Unique[1],Unique[x])}))/Range
  strand <- unlist(apply(matrix(as.character(Unique$SYMBOL),ncol=1),1,function(x){subset <- Gene.range[Gene.range$SYMBOL %in% x]
                                                                                  return(unique(as.character(strand(subset)))[1])}))
  
  df <- data.frame(Symbols=as.character(Unique$SYMBOL),
                   x1=Distance,
                   strand = strand,
                   colors=rep("black",length(Unique)),
                   label= rep("No",length(Unique)),
                   stringsAsFactors = F)
                   
  df$colors[df$Symbols %in% special$Symbols] <- special$colors
  df$label[df$Symbols %in% special$Symbols] <- "Yes"
  if(save) pdf(paste0(fn,".pdf"))
  vp <- viewport(h=5, w=6)
  grid.lines(c(0.2,1),c(0.3,0.3),gp=gpar(lwd=3))
  grid.text(fn,0.1,0.3,gp=gpar(lwd=3))
  grid.rect(targetD,0.3,gp=gpar(col=target.col,fill=target.col),width=0.015,height=0.01)
  .plot.arrow <- function(x){
    start <- as.numeric(x["x1"])
    if(x["strand"] == "-"){
      add <- -0.025
      if(x["label"] == "Yes"){
        grid.curve(x1=start, y1=0.3, x2=start+add*2, y2=0.35, arrow=arrow(length = unit(0.02,"npc")),shape=0,gp=gpar(col=x["colors"],lwd=3))
        grid.lines(x=c(start,start+0.05),
                   y=c(0.36,0.5),gp=gpar(lwd=2))
        grid.text(x["Symbols"],start+0.05,0.52,gp=gpar(lwd=3))
      }else{
        grid.curve(x1=start, y1=0.3, x2=start+add, y2=0.325, arrow=arrow(length = unit(0.01,"npc")),shape=0,gp=gpar(col=x["colors"],lwd=3))
      }
    }else{
      add <- 0.025
      if(x["label"] == "Yes"){
        grid.curve(x1=start, y1=0.3, x2=start+add*2, y2=0.35, arrow=arrow(length = unit(0.02,"npc")),shape=0,gp=gpar(col=x["colors"],lwd=3),curvature=-1)
        grid.lines(x=c(start,start+0.05),
                   y=c(0.36,0.5),gp=gpar(lwd=2))
        grid.text(x["Symbols"],start+0.05,0.52,gp=gpar(lwd=3))
      }else{
        grid.curve(x1=start, y1=0.3, x2=start+add, y2=0.325, arrow=arrow(length = unit(0.01,"npc")),shape=0,gp=gpar(col=x["colors"],lwd=3),curvature=-1)
      }
    }
   
  }
  for(i in 1:nrow(df)){
    .plot.arrow(as.matrix(df)[i,])
  }
  if(save) dev.off()
}

