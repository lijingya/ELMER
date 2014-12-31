## common colors
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
redGreen <- colorRampPalette(c("green","black","red"))

#'cluster functions
#'@param x A matrix 
#'@param Rowv A boolean determines if the row dendrogram should be computed and reordered. 
#'@param Colv A boolean determines if the column dendrogram should be computed and reordered.
#'@param distfun function used to compute the distance (dissimilarity) between both rows and columns. Defaults to dist.
#'@param distMethod A character to specify method for computing distance. Default to euclidean. See detail for other methods.
#'@param hclustfun function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms. Defaults to fastcluster::hclust. Should take as argument a result of distfun and return an object to which as.dendrogram can be applied.
#'@param hclustMethod A character to specify method for computing clustering. Default to complete
#'@param Distance.row A vector of distance value for rows. If Distance.row was specified, distance calculation step will be skiped for rows.
#'@param Distance.col A vector of distance value for columns. If Distance.col was specified, distance calculation step will be skiped for columns.
#'@return A list contains: x the original matrix; rowInd order of row after clustering; ddr dendrograms for rows; colInd order of columns after clustering; ddc dendrograms for columns.
#'@details distMethod euclidean: Usual square distance between the two vectors (2 norm).maximum: Maximum distance between two components of x and y (supremum norm).manhattan: Absolute distance between the two vectors (1 norm). canberra: sum(|x_i - y_i| / |x_i + y_i|). Terms with zero numerator and denominator are omitted from the sum and treated as if the values were missing. This is intended for non-negative values (e.g. counts): taking the absolute value of the denominator is a 1998 R modification to avoid negative distances.binary: (aka asymmetric binary): The vectors are regarded as binary bits, so non-zero elements are ‘on’ and zero elements are ‘off’. The distance is the proportion of bits in which only one is on amongst those in which at least one is on. minkowski: The p norm, the pth root of the sum of the pth powers of the differences of the components.
#'@details hclustMethod: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
cluster.main <- function(x, Rowv = TRUE, Colv = TRUE, distfun = dist, distMethod="euclidean", hclustfun = fastcluster::hclust, hclustMethod="complete", 
                         Distance.row=NULL, Distance.col=NULL){
  ##use faster cluster package
  library(fastcluster)
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("'x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("'x' must have at least 2 rows and 2 columns")
  if (Rowv){
    if(is.null(Distance.row)){
      Distance.row <- distfun(x,method=distMethod)
    }
    hcr <- hclustfun(Distance.row,method=hclustMethod)
    ddr <- as.dendrogram(hcr)
    if (nr != length(rowInd <- order.dendrogram(ddr))){
      stop("row dendrogram ordering gave index of wrong length")
    }  
  }else{
    rowInd <- 1:nr
    ddr <- NULL
  }
  
  if(Colv){
    if(is.null(Distance.col)){
      Distance.col <- distfun(t(x),method=distMethod)
    }
    hcc <- hclustfun(Distance.col,method=hclustMethod)
    ddc <- as.dendrogram(hcc)
    if (nc != length(colInd <- order.dendrogram(ddc))){
      stop("col dendrogram ordering gave index of wrong length")
    }  
  }else{
    colInd <- 1:nc
    ddc <- NULL
  }
  
#  x <- x[rowInd, colInd]
  out <- list(x=x,rowInd = rowInd, ddr=ddr,colInd = colInd,ddc=ddc,hcc=hcc)
  return(out)
}
 

#' output heatmap
#' @param x output from cluster.main.
#' @param margins a character vector of variable names to compute margins for. 
#' @param labRow a character vector of labels for rows of matrix in x.
#' @param labCol a chracter vector of lables for columns of matrix in x.
#' @param nonlab a boolean to determine no labels for rows and column.
#' @param nonlab.row a boolean to determine no labels for rows.
#' @param nonlab.row a boolean to determine no labels for columns.
#' @param ... parameters for image function.
#' @return A heatmap
heatmap.main <- function(x,margins=c(5,0.5,0.5,5), labRow=NULL, labCol=NULL, nonlab=F,nonlab.row=F,nonlab.col=F, xlab=NULL, ylab=NULL, col=heat.colors(225), zlim=NULL, cexRow = 0.2 +1/log10(nr), cexCol = 0.2 + 1/log10(nc)){
  if (length(di <- dim(x$x)) != 2 || !is.numeric(x$x)) stop("'x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) stop("'x' must have at least 2 rows and 2 columns")
#labcol, labrow,  
  if (is.null(labRow)) {
    if (is.null(rownames(x$x))) {
      labRow <- (1:nr)[x$rowInd]
    }else{
      labRow <-rownames(x$x)[x$rowInd]
    }    
  }else{
    labRow <- labRow[x$rowInd]
  } 
  if (is.null(labCol)){
    if (is.null(colnames(x$x))){
      labCol <- (1:nc)[x$colInd]
    }else{
      labCol <-colnames(x$x)[x$colInd]
    }  
  }else{
    labCol <- labCol[x$colInd]
  } 
  if(nonlab){
    labRow<- NULL
    labCol <- NULL
  }
  x$x <- x$x[x$rowInd, x$colInd]
  par(mar=margins)
  if (is.null(zlim)){
    image(1:nc, 1:nr, t(x$x), xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col)
  }else{
    image(1:nc, 1:nr, t(x$x), zlim=zlim, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col)
  }
  
  if(!(nonlab.col|nonlab)){
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
  }
  
  if (!is.null(xlab)){
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  }  
  
  if (!(nonlab.row|nonlab)){
    axis(2, 1:nr, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
  }
 
  
  if (!is.null(ylab)){
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  }     
}

#' making side bars for heatmap
#' @param x A matrix which is same order of the original matrix before cluster. See details
#' @param side A character which are either 'colside' or 'rowside' specifying where the side bar locates.
#' @param order A vector of number specifying order of side bars. 
#' @param ... parameters for image function.
#' @return a side bar for heatmap
#' @details x must be the matrix. If it is colbars, the row number of col bars should be the same as the col number of the matrix. If it is rowbars, the row number of row bars should be the same as the row number of the matrix
side.bars <- function(x,side="colside",order=NULL,margins=c(5,0.5,0.5,5),lab=NULL,col=heat.colors(225),zlim=NULL,cexlab = 0.2 +1/log10(nr)){
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) stop("'x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  names <- colnames(x)
  par(mar=margins)
  if(side == "colside"){
    cbar <- x[order,]
    cbar=matrix(cbar, nrow=nr)
    if (is.null(zlim)){
      image(cbar, col = col, axes = FALSE) 
    }else{
      image(cbar, col = col, zlim=zlim, axes = FALSE) 
    }
    
    if (is.null(lab)) {
      if(nc==1){
        axis(2, 0 , colnames(cbar), las=2, tick=FALSE,cex.axis=cexlab)
      }else{
        axis(2, 0:(dim(cbar)[2]-1)/(dim(cbar)[2]-1) , names, las=2, tick=FALSE,cex.axis=cexlab)
      }
      
    }else{
      if(nc==1){
        axis(2, 0, lab, las=2, tick=FALSE,cex.axis=cexlab)
      }else{
        axis(2, 0:(dim(cbar)[2]-1)/ (dim(cbar)[2]-1), lab, las=2, tick=FALSE,cex.axis=cexlab)
      }   
    }
  }
  if(side == "rowside"){
    rbar <- x[order,]
    rbar<- matrix(rbar, nrow=nr)
    if (is.null(zlim)){
      image(t(rbar), col = col, axes = FALSE) 
    }else{
      image(t(rbar), col = col, zlim=zlim, axes = FALSE) 
    }
    
    if (is.null(lab)) {
      if(nc==1){
        axis(1, 0 , colnames(rbar), las=2, tick=FALSE,cex.axis=cexlab)
      }else{
        axis(1, 0:(dim(rbar)[2]-1)/ (dim(rbar)[2]-1) , names, las=2, tick=FALSE,cex.axis=cexlab)
      }
      
    }else{
      if(nc==1){
        axis(1, 0, lab, las=2, tick=FALSE,cex.axis=cexlab)
      }else{
        axis(1, 0:(dim(rbar)[2]-1)/ (dim(rbar)[2]-1), lab, las=2, tick=FALSE,cex.axis=cexlab)
      }   
    }  
  }
  box("plot",col="black")
  #add boader
  if(nc>1){
    
    for(i in 1:nc-1){
      if(side=="colside"){
        abline(h=par("usr")[3]+i*(par("usr")[4]-par("usr")[3])/nc,col="black")
      }else{
        abline(v=par("usr")[1]+i*(par("usr")[2]-par("usr")[1])/nc,col="black")
      }
      
    }
  }
}


#' making side bars for heatmap using segment to solve the problem in image that line will very tiny when a lot of samples.
#' @param x A matrix which is same order of the original matrix before cluster. See details
#' @param side A character which are either 'colside' or 'rowside' specifying where the side bar locates.
#' @param order A vector of number specifying order of side bars. 
#' @param ... parameters for image function.
#' @return a side bar for heatmap
#' @details x must be the matrix. If it is colbars, the row number of col bars should be the same as the col number of the matrix. If it is rowbars, the row number of row bars should be the same as the row number of the matrix
side.bars2 <- function(x,side="colside",order=NULL,margins=c(5,0.5,0.5,5),lab=NULL,zlim=NULL,cexlab = 0.2 +1/log10(nr)){
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) stop("'x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  names <- colnames(x)
  par(mar=margins)
  if(side == "colside"){
    cbar <- x[order,]
    cbar=matrix(cbar, nrow=nr)
    image(cbar, col = "white", axes = FALSE)
    dis <- par("usr")[4]-par("usr")[3]
    box("plot",col="black")
    #add boader
    if(nc>1){
      for(i in 1:nc){
        abline(h=par("usr")[3]+i*dis/nc,col="black")
        for(nn in which(cbar[,i]==1)){
          segments(nn/nr,par("usr")[3]+(i-1)*dis/nc,nn/nr,par("usr")[3]+i*dis/nc)
        }
      }
    }else{
      for(nn in which(cbar==1)){
        segments(nn/nr,par("usr")[3],nn/nr,par("usr")[3]+dis/nc)
      }
    }
    if (!is.null(lab)) {
     names <- labs 
    }
    if(nc==1){
      axis(2, 0 , colnames(cbar), las=2, tick=FALSE,cex.axis=cexlab)
    }else{
      axis(2, 0:(dim(cbar)[2]-1)/(dim(cbar)[2]-1) , names, las=2, tick=FALSE,cex.axis=cexlab)
    }
  }
  if(side == "rowside"){
    rbar <- x[order,]
    rbar<- matrix(rbar, nrow=nr)
    image(t(rbar), col = "white", axes = FALSE)
    dis <- par("usr")[2]-par("usr")[1]
    box("plot",col="black")
    #add boader
    if(nc>1){
      for(i in 1:nc){
        abline(v=par("usr")[1]+i*dis/nc,col="black")
        for(nn in which(rbar[,i]==1)){
          segments(par("usr")[1]+(i-1)*dis/nc,nn/nr,par("usr")[1]+i*dis/nc,nn/nr,)
        }
      }
    }else{
      for(nn in which(rbar==1)){
        segments(par("usr")[1],nn/nr,par("usr")[1]+dis/nc,nn/nr)
      }
    }
    if (!is.null(lab)) {
      names <- labs 
    }  
    if(nc==1){
      axis(1, 0 , colnames(rbar), las=2, tick=FALSE,cex.axis=cexlab)
    }else{
      axis(1, 0:(dim(rbar)[2]-1)/ (dim(rbar)[2]-1) , names, las=2, tick=FALSE,cex.axis=cexlab)
    }
  }
  
}



#' keyplot
#' @param x A matrix.
#' @param col A vector of colors to define colors.
#' @param breaks a vector of number to define how many label in x axis.
#' @param extremes a vector of number to define the edge value of the key. Default is NULL.
#' @param ... parameters for image.
#' @return A keyplot.
keyplot <- function(x,col,breaks,extremes=NULL,texts,margin=c(1,1,1,1),cex.axis=par("cex.axis")){
  #x axis lable position: axis(at=..)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  
  #break is to set up the x label .
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  
  if (is.null(extremes)){ 
    breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),length = breaks) 
    #get the minimul and maximal data
    min.raw <- min(x, na.rm = TRUE)
    max.raw <- max(x, na.rm = TRUE)
  }else{
    breaks <- seq(extremes[1], extremes[2], length = breaks)
    min.raw <- extremes[1]
    max.raw <- extremes[2]
  }
  par(mar=margin)
  # color sequential data
  z <- seq(min.raw, max.raw, length = length(col))
  image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
        xaxt = "n", yaxt = "n")
  #add the x axis label
  lv <- pretty(breaks)
  xv <- scale01(as.numeric(lv), min.raw, max.raw)
  axis(1, at = xv, labels = lv,cex.axis=cex.axis)
  
  #add label
  if(!missing(texts)){
    mtext(side=1,texts,line=2)
  }else{
    mtext(side=1,"Value",line=2)
  }
  title("Color Key",cex.main=2)  
}


#' Plot the dendro tree.
#' @param Rdend dendrograms for rows.
#' @param Cdend dendrograms for columns.
#' @param title The main title (on top).
#' @param cex.title A numerical value giving the amount by which plotting title text and symbols should be magnified relative to the default. This starts as 1 when a device is opened, and is reset when the layout is changed.
#' @return a graph of dendrograms tree.
dendro.plot <- function(Rdend=NULL, Cdend=NULL,margin=c(0.5,0.5,0.5,0.5),title=NULL,cex.title){
  par(mar = margin)
  if (!is.null(Rdend)) {
    plot(Rdend, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  if (!is.null(Cdend)){
    plot(Cdend, axes = FALSE, xaxs = "i", leaflab = "none")
  } 
  if(!is.null(title)){
    title(main=title,cex=cex.title)
  }
}

mat=rbind(c(NA,NA,7),c(6,NA,5),c(NA,NA,4),c(3,2,1))


#add more information to the figure
heat.info <- function(main=NULL,cex=1.5,mars=c(4,4,4,4)){
  plot.new()
  par(mar=mars)
  text(par("usr")[2],par("usr")[4], main, cex = cex)
}

#' add lengend
#' @param Labels A vector of characters
#' @param cols A vector of colors for each characters in Labels
#' @param ... parameters in image function.
#' @return A legend
AddLegend <- function(Labels=NULL,cols=NULL,margins=c(1,1,1,1),lab.las=2,cexlab=0.2 +1/log10(length(cols))){
  par(mar=margins)
  if(is.null(Labels)) Labels <- 1:length(cols)
  z <- 1:length(cols)
  image(matrix(z,ncol=1), col = cols, axes = FALSE)
  box("plot",col="black")
  labposition <- c()
  for(i in 1:length(cols)-1){
    abline(v=par("usr")[1]+i*(par("usr")[2]-par("usr")[1])/length(cols),col="black") 
    labposition <- c(labposition,i*(par("usr")[2]-par("usr")[1])/length(cols))
  }
  axis(1, labposition, Labels, las=lab.las, tick=FALSE,cex.axis=cexlab)  
}


#--------------color function-------------------------------
#brewer.pal(n, name) :makes the color palettes from ColorBrewer available as R palettes.
#display.brewer.pal(n, name):  displays the selected palette in a graphics window.
#display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE): displays the a few palettes simultanueously in a graphics window.
#There are 3 types of palettes, sequential, diverging, and qualitative.
#1. Sequential palettes are suited to ordered data that progress from low to high. Lightness steps dominate the look of these schemes, with light colors for low data values to dark colors for high data values. 
#2. Diverging palettes put equal emphasis on mid-range critical values and extremes at both ends of the data range. The critical class or break in the middle of the legend is emphasized with light colors and low and high extremes are emphasized with dark colors that have contrasting hues. 
#3. Qualitative palettes do not imply magnitude differences between legend classes, and hues are used to create the primary visual differences between classes. Qualitative schemes are best suited to representing nominal or categorical data.

#The sequential palettes names are 
#Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd

#All the sequential palettes are available in variations from 3 different values up to 9 different values.

#The diverging palettes are 
#BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral

#All the diverging palettes are available in variations from 3 different values up to 11 different values.

#For qualitative palettes, the lowest number of distinct values available always is 3, but the largest number is different for different palettes. It is given together with the palette names in the following table.

#Accent   8
#Dark2   8
#Paired	 12
#Pastel1	 9
#Pastel2	 8
#Set1	 9
#Set2	 8
#Set3	 12

# GenerateColor <- function(n,name){
#   if (!require("RColorBrewer")) {
#     install.packages("RColorBrewer")
#     library(RColorBrewer)
#   }
#   cols <- colorRampPalette(brewer.pal(n,name))
#   return(cols)
# }
# 
# Color2Num <- function(x){
#   num <- unique(x)
#   Numbers <- rep(0,length(x))
#   count <- 0
#   for(i in num){
#     Numbers[x %in% i] <-count
#     count <- count+1
#   }
#   return(Numbers)
# }
# 
# 
# #-------------convert number to color -----------------------
# #x is vector
# SetcolorNum <- function (x,TotalNum=255, nameBrewer,Add=0,extremes=NULL){
#   cols<- GenerateColor(9,nameBrewer)
#   col <- cols(TotalNum)
#   breaks <- length(col) + 1
#   if (is.null(extremes)) 
#     breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),length = breaks)
#   else breaks <- seq(extremes[1], extremes[2], length = breaks)
#   Cuts <- cut(x,breaks=breaks,labels=FALSE)
#   Cuts <- Add+Cuts
#   return(Cuts) 
# }


#' Normalization to 0 to 1
#' @param x A matrix.
#' @param col A boolean to determine normalize by column or not.
#' @param row A boolean to determine normalize by row or not.
#' @param na.rm A boolean to determine to remove na number or not.
#' @return A normalized matrix.
Normalize <- function (x,col=FALSE,row=FALSE,na.rm=FALSE){
  if(col){
    ColMax <- apply(x,2,max,na.rm=na.rm)
    ColMin <- apply(x,2,min,na.rm=na.rm)
    x <- t((t(x)-ColMin)/ColMax)
  }
  if(row){
    RowMax <- apply(x,1,max,na.rm=na.rm)
    RowMin <- apply(x,1,min,na.rm=na.rm)
    x <- (x-RowMin)/RowMax
  }
  
  return(x)
}

#' Normalization based on mean
#' @param x A matrix.
#' @param col A boolean to determine normalize by column or not.
#' @param row A boolean to determine normalize by row or not.
#' @param na.rm A boolean to determine to remove na number or not.
#' @return A normalized matrix.
NormalizeMean <- function (x,col=FALSE,row=FALSE,na.rm=FALSE){
  
  if(col){
    Mean <- colMeans(x,na.rm=na.rm)
    SD <- apply(x,2,sd,na.rm=na.rm)
    x <- t((t(x)-Mean)/SD)
    x[,SD==0] <- 0 
  }
  if(row){
    Mean <- rowMeans(x,na.rm=na.rm)
    SD <- apply(x,1,sd,na.rm=na.rm)
    x <- (x-Mean)/SD
    x[SD==0,] <- 0
  }
  
  return(x)
}

#' Normalization based on median
#' @param x A matrix.
#' @param col A boolean to determine normalize by column or not.
#' @param row A boolean to determine normalize by row or not.
#' @param na.rm A boolean to determine to remove na number or not.
#' @return A normalized matrix.
NormalizeMedian <- function (x,col=FALSE,row=FALSE,na.rm=FALSE){
  if(col){
    Median <- apply(x,2,median,na.rm=na.rm)
    x <- t((t(x)-Median))
  }
  if(row){
    Median <- apply(x,1,median,na.rm=na.rm)  
    x <- x-Median
  }
  
  return(x)
}

#' binary data
#' @param x A matrix.
#' @param Break A value to binarize the data.
#' @param Break2 A value to cut value to 3 categories.
#' @return A binarized matrix.
Binary <- function(x,Break=0.3,Break2=NULL){
  if(!is.numeric(x)) stop("x need to be numeric") 
  change <- x
  if(is.null(Break2)){
    change[x > Break] <- 1
    change[x < Break | x== Break] <- 0
  }else{
    change[x < Break | x== Break] <- 0
    change[x> Break & x < Break2] <- NA
    change[x > Break2 | x== Break2] <-1 
  }
  
  return(change)    
}

# 
# ##make multisidebars-------------------------------------------------------
# MultiSide.bars <- function(data=list(),side="colside",order=NULL,margins=c(5,0.5,0.5,5),col=NULL,zlim=NULL,cexlab = 0.2 +1/log10(nn)){
#   nn <- length(data)
#   if (is.null(order)) stop ("order should not be NULL")
#   if(is.null(zlim)) zlim <- list()
#   if(is.null(col)){
#     Col <- NULL
#     col <- list()
#   }else{
#     Col <- TRUE
#   }
#   for (i in names(data)){
#     print(i)
#     if(is.null(zlim)) zlim[[i]]<- NULL 
#     Bar.name <- i
#     One <- data[[i]]
#     if(is.numeric(One)){
#       if(is.null(Col)) col[[i]] <- jet.colors(255)
#       side.bars(matrix(One,ncol=1),side=side,order=order,margins=margins,col=col[[i]],lab=Bar.name,cexlab=cexlab,zlim=zlim[[i]])
#       keyplot (One,col=col[[i]],margin=c(1,1,0.5,1),extremes=zlim[[i]])
#     }else if(is.factor(One)){
#       if(is.null(Col)) col[[i]] <- 1:length(unique(One))
#       side.bars(matrix(as.numeric(One),ncol=1),side=side,order=order,margins=margins,col=col[[i]][sort(unique(as.numeric(One)))],lab=Bar.name,cexlab=cexlab)
#       AddLengend (levels(One)[sort(unique(as.numeric(One)))],cols=col[[i]][sort(unique(as.numeric(One)))],margins=c(1,1,0.5,1),cexlab=cexlab,lab.las=1)
#     }else if(is.character(One)){
#       if(is.null(Col)) col[[i]] <- 1:length(unique(One))
#       One <- factor(One)
#       side.bars(matrix(as.numeric(One),ncol=1),side=side,order=order,margins=margins,col=col[[i]][sort(unique(as.numeric(One)))],lab=Bar.name,cexlab=cexlab)
#       AddLengend (levels(One)[sort(unique(as.numeric(One)))],cols=col[[i]][sort(unique(as.numeric(One)))],margins=c(1,1,0.5,1),cexlab=cexlab,lab.las=1)
#     }
#   }
# }



#' lable linear regression formula 
#' @param df A data.frame object contains two variables: dependent variable (Dep) and explanation variable (Exp).
#' @return a linear regression formula
lm_eqn = function(df){
  m = lm(Dep ~ Exp, df);
  eq <- substitute(italic(y) == a + (b) %.% italic(x)*"\n"~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}



##PeakToVenn
#Peaks :List the GRange format objects.
#... : the findOverlap options.
#the overlap parts belongs to the peaks set before.
#' Making peak sets venn diagram
#' @param Peaks A list of Peak sets.
#' @param ... parameters from VennDiagram package.
#' @return A venn diagram of peaks.
PeakToVenn <- function(Peaks,...){
  library(VennDiagram)
  if(length(Peaks)==2){
    P1 <- Peaks[[1]]
    P2 <- Peaks[[2]]
    values(P1) <- NULL
    values(P2) <- NULL
    Over <- findOverlaps(P1,P2)
    A <- 1: length(P1)
    OverNum <- length(unique(queryHits(Over)))
    B_unique <- length(P2)-length(unique(subjectHits(Over)))
    B <- (length(P1)-OverNum+1):(length(P1)-OverNum+B_unique)
    tmp <- list(A,B)
    names(tmp) <- names(Peaks)
    venn.plot <- venn.diagram(
      x = tmp,
      filename = NULL,
      fill = c("cornflowerblue", "darkorchid1"),
      alpha = 0.75,
      label.col = "black",
      fontfamily = "serif",
      fontface = "bold",
      cat.col = c("cornflowerblue", "darkorchid1"),
      cat.fontfamily = "serif",
      cat.fontface = "bold",
      cat.dist = c(0.03, 0.03),
      cat.pos = c(-20, 14),
      ...
    );
  }else if (length(Peaks)==3){
    for (i in 1:length(Peaks)){
      values(Peaks[[i]]) <- NULL
    }
    ## n1: 1 unique, n2: 2 unique, n3: 3 unique
    ##n123_1: in n123 what is the number of 1. n123_2 ...of 2, n123_3 ...of 3
    ##n12_1: in n12 what is the number of 1.  n12_2 ...of 2
    ##n13_1: in n13 what is the number of 1.  n13_3 ...of 3
    ##n23_2: in n23 what is the number of 2.  n23_3 ...of 3
    ##n1
    over1 <- findOverlaps(Peaks[[1]],c(Peaks[[2]],Peaks[[3]]))
    n1 <- length(Peaks[[1]]) - length(unique(queryHits(over1)))
    ##n2
    over1 <- findOverlaps(Peaks[[2]],c(Peaks[[1]],Peaks[[3]]))
    n2 <- length(Peaks[[2]])-length(unique(queryHits(over1)))
    ##n3
    over1 <- findOverlaps(Peaks[[3]],c(Peaks[[2]],Peaks[[1]]))
    n3 <- length(Peaks[[3]])-length(unique(queryHits(over1)))
    ##n123
    tmpPeak <- Peaks[[1]][unique(queryHits(findOverlaps(Peaks[[1]],Peaks[[2]])))]
    n123_1 <- length(unique(queryHits(findOverlaps(tmpPeak,Peaks[[3]]))))
    tmpPeak <- Peaks[[2]][unique(queryHits(findOverlaps(Peaks[[2]],Peaks[[1]])))]
    n123_2 <- length(unique(queryHits(findOverlaps(tmpPeak,Peaks[[3]]))))
    tmpPeak <- Peaks[[3]][unique(queryHits(findOverlaps(Peaks[[3]],Peaks[[1]])))]
    n123_3 <- length(unique(queryHits(findOverlaps(tmpPeak,Peaks[[2]]))))
    #n12
    n12_1 <- length(unique(queryHits(findOverlaps(Peaks[[1]],Peaks[[2]]))))-n123_1
    n12_2 <- length(unique(queryHits(findOverlaps(Peaks[[2]],Peaks[[1]]))))-n123_2
    #n13
    n13_1 <- length(unique(queryHits(findOverlaps(Peaks[[1]],Peaks[[3]]))))-n123_1
    n13_3 <- length(unique(queryHits(findOverlaps(Peaks[[3]],Peaks[[1]]))))-n123_3
    #n23
    n23_2 <- length(unique(queryHits(findOverlaps(Peaks[[2]],Peaks[[3]]))))-n123_2
    n23_3 <- length(unique(queryHits(findOverlaps(Peaks[[3]],Peaks[[2]]))))-n123_3
    tmp <- list()
    tmp[[1]] <- 1:length(Peaks[[1]])
    tmp[[2]] <- c(1:n123_1,              ##common sets
                  (n123_1+1):(n123_1+n12_1),  #n12_1
                  (length(Peaks[[1]])+1):(length(Peaks[[1]])+n2+n23_2))  ##unique to 2
    tmp[[3]] <- c(1:n123_1,
                  (n123_1+n12_1+1):(n123_1+n12_1+n13_1),      ##n13_1
                  (length(Peaks[[1]])+1):(length(Peaks[[1]])+n23_2),      #n23_2
                  (max(tmp[[2]])+1):(max(tmp[[2]])+n3))
    names(tmp) <- names(Peaks)

    subtile <- sprintf("1:%s,2:%s,3:%s\nn1:%d,n2:%d,n3%d,n12_1:%d,n12_2:%d,n13_1:%d,\nn13_3:%d,n23_2:%d,n23_3:%d,n123_1:%d,n123_2:%d,n123_3:%d",names(tmp)[1],names(tmp[2]),names(tmp)[3],
                       n1,n2,n3,n12_1, n12_2,n13_1,n13_3,n23_2,n23_3,n123_1,n123_2,n123_3)
    venn.plot <- venn.diagram(
      x = tmp,
      filename = NULL,
      col = "transparent",
      fill = c("red", "blue", "green"),
      alpha = 0.5,
      label.col = c("darkred", "white", "darkblue", "white",
                    "white", "white", "darkgreen"),
      fontfamily = "serif",
      fontface = "bold",
      cat.default.pos = "text",
      cat.col = c("darkred", "darkblue", "darkgreen"),
      cat.fontfamily = "serif",
      cat.dist = c(0.06, 0.06, 0.03),
      cat.pos = 0,
      sub=subtile,
      ...
    );
  }
  return (venn.plot)
}

# convert chr matrix to color ---------------------------------------------
##x is matrix with charatcter.

ChrToColor <- function(x,cols=NULL,simple=FALSE){
  if(class(x) %in% "data.frame") x <- as.matrix(x)
  Levels <- sort(unique(as.vector(x)))
  if(any(is.na(as.vector(x)))) Levels <- c("missing",Levels)
  if(is.null(cols)) cols <- 1:length(Levels)
  out <- mat.or.vec(nr=nrow(x),nc=ncol(x))
  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)
  if(simple){
    for(i in 1:length(Levels)){
      if(Levels[i] %in% "missing") out[is.na(x)] <- cols[i]
      out[x %in% Levels[i]] <- cols[i]
    }
    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    return(out)
  }else{
    for(i in 1:length(Levels)){
      if(Levels[i] %in% "missing") out[is.na(x)] <- cols[i]
      out[x %in% Levels[i]] <- i
    }
    colnames(out) <- colnames(x)
    rownames(out) <- rownames(x)
    tmp <- list(x=out,cols=cols,levels=Levels)
    return(tmp)
  }
}



# ## smooth DNA methylation matrix
# Smooth <- function(x, bin=10, cores=6){
#   nr <- dim(x)[1]
#   nc <- dim(x)[2]
#   start <- 1
#   end <- nc-bin   
#   cl <- makeCluster(cores,"SOCK")
#   out <- parSapplyLB(cl, start:end, function(x, data.matrix, bin){sub.matrix <- data.matrix[,x:(x+bin-1)]
#                                                                   out <- rowMeans(sub.matrix,na.rm=T)
#                                                                   return(out)}, data.matrix=x, bin=bin, simplify=F)
#   stopCluster(cl)
#   out <- do.call(cbind,out)
#   return(out)
# }
# 
# # bin DNA methylation matrix
# Bin <- function(x, size, cores=6){
#   nr <- dim(x)[1]
#   nc <- dim(x)[2]
#   isInt <- function(n) {
#     return (ceiling(n) == n);
#  }
#  if(!isInt(size/2)) size <- 1 + size
#  if(isInt(nc/2)){
#    center = nc/2+1
#    k = floor((nc-center-size/2+1)/size)
#    Breakpoints <- as.vector(unlist(apply(matrix(0:k,ncol=1), 1,function(x,center,size){ upstream <- center+size/2-1+x*size
#                                                          downstream <- center - size/2 - x*size
#                                                          out <- c(upstream, downstream)
#                                                          return(out)}, center=center, size=size)))
#    if(nc- center - size/2 + 1 - k*size > size/2) Breakpoints <- c(1,nc)
#    Breakpoints <- sort(Breakpoints)
#    Sites <- Breakpoints
#    Sites[Sites < center] <- Sites[Sites < center] - center+1
#    Sites[Sites > center] <- Sites[Sites > center] - center
#  }else{
#    center = ceiling(nc/2)
#    k = floor((nc-center-size/2)/size)
#    Breakpoints <- as.vector(unlist(apply(matrix(0:k,ncol=1), 1,function(x,center,size){ upstream <- center+size/2+x*size
#                                                                                         downstream <- center - size/2 - x*size
#                                                                                         out <- c(upstream, downstream)
#                                                                                         return(out)}, center=center, size=size)))
#    if(nc- center - size/2 - k*size > size/2) Breakpoints <- c(1,nc)
#    Breakpoints <- sort(Breakpoints)
#    Sites <- Breakpoints-center
#  }
#  Sites <- sapply(1:(length(Sites)-1), function(x){ out <- mean(Sites[c(x,x+1)])
#                                                    return(out)})
#  cl <- makeCluster(cores,"SOCK")
#  out <- parSapplyLB(cl, 1:(length(Breakpoints)-1), function(x, data.matrix, Breakpoints){sub.matrix <- data.matrix[,Breakpoints[x]:(Breakpoints[x+1]-1)]
#                                                                                           out <- rowMeans(sub.matrix,na.rm=T)
#                                                                                           return(out)}, data.matrix=x, Breakpoints=Breakpoints, simplify=F)
#  stopCluster(cl)
#  out <- list(matrix=do.call(cbind,out), sites=Sites)
#  return(out)
# }