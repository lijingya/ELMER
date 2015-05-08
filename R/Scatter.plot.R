#'scatter.plot
#'@importFrom ggplot2 ggsave
#'@param mee A MEE.data object includes DNA methylation data, expression data, 
#'probeInfo and geneInfo.
#'@param byPair A list: byPair =list(probe=c(),gene=c()); probe contains a vector 
#'of probes' name and gene contains a vector of gene ID. The length of probe 
#'should be the same with length of gene. Output see detail.
#'@param byProbe A list byProbe =list(probe=c(), geneNum=20); probe contains 
#'a vector of probes'name and geneNum specify the number of gene near the probes 
#'will ploted. 20 is default for geneNum. Output see detail.
#'@param byTF A list byTF =list(TF=c(), probe=c()); TF contains a vector of TF's 
#'symbol and probe contains the a vector of probes' name. Output see detail.
#'@param category A vector labels subtype of samples or a character which is the 
#'column name in the sampleInfo in the MEE.data object. Once specified, samples 
#'will label different color. The color can be customized by using color.value. 
#'@param dir.out A path specify the directory to which the figures will be saved. 
#'Current directory is default.
#'@param save A logic. If true, figure will be saved to dir.out.
#'@param ... color.value, lm_line in scatter function
#'@details byPair The output will be scatter plot for individual pairs.
#'@details byProbe The output will be scatter plot for the probe and nearby genes.
#'@details byTF The output will be scatter plot for the TFs and the average 
#'DNA methylation at the probes set specified in byTF list.
#'@return Scatter plots.
#'@export
#'@examples
#'load(system.file("extdata","mee.example.rda",package = "ELMER"))
#'scatter.plot(mee,byProbe=list(probe=c("cg19403323"),geneNum=20), 
#'            category="TN", save=FALSE)
#'scatter.plot(mee,byProbe=list(probe=c("cg19403323"),geneNum=20), 
#'            category="TN", save=TRUE) ## save to pdf
#'# b. generate one probe-gene pair
#'scatter.plot(mee,byPair=list(probe=c("cg19403323"),gene=c("ID255928")),
#'             category="TN", save=FALSE,lm_line=TRUE) 
            
            
scatter.plot <- function(mee,byPair=list(probe=c(),gene=c()),
                         byProbe=list(probe=c(),geneNum=20),byTF=list(TF=c(),probe=c()), 
                         category=NULL, dir.out ="./", save=TRUE, ...){
  if(missing(mee)) stop("A MEE.data object should be included.")
  if(!is.null(category) && length(category)==1) 
    category <- getSample(mee,cols = category)
  if(length(byPair$probe) != 0){
    if(length(byPair$probe)!=length(byPair$gene))
      stop("In pairs, the length of probes should be the same with the length of genes.")
    for(i in 1:length(byPair$probe)){
      probe <- byPair$probe[i]
      gene <- byPair$gene[i]
      symbol <- getSymbol(mee,geneID=gene)
      P <- scatter(meth=getMeth(mee,probe=probe),
                   exp=getExp(mee,geneID = gene ),category=category, 
              xlab=sprintf("DNA methyation at %s",probe), 
              ylab=sprintf("%s gene expression",symbol), 
              title=sprintf("%s_%s",probe,symbol),
              ...)
      if(save) ggsave(filename = sprintf("%s/%s_%s.bypair.pdf",dir.out,probe,symbol),
                      plot = P,useDingbats=FALSE, width=7, height = 6)
    }
  }
  
  if(length(byProbe$probe)!=0){
    nearGenes <- GetNearGenes(TRange=getProbeInfo(mee,probe=byProbe$probe),
                              geneAnnot=getGeneInfo(mee),geneNum = byProbe$geneNum)
    for(i in byProbe$probe){
      probe <- i
      gene <- nearGenes[[i]]$GeneID
      symbol <- getSymbol(mee,geneID=gene)
      exp <- getExp(mee,geneID=gene)
      rownames(exp) <- symbol
      P <- scatter(meth=getMeth(mee,probe=probe), exp=exp,category=category,
                   xlab=sprintf("DNA methyation at %s",probe), 
                   ylab=sprintf("Gene expression"), 
                   title=sprintf("%s nearby %s genes",probe,byProbe$geneNum),
                   ...)
      if(save) ggsave(filename = sprintf("%s/%s.byprobe.pdf",dir.out,probe),
                      plot = P,useDingbats=FALSE)
    }
  }
  
  if(length(byTF$TF)!=0){
    meth <- colMeans(getMeth(mee,probe=byTF$probe),na.rm = TRUE)
    gene <- getGeneID(mee,symbol=byTF$TF)
    exp <- getExp(mee,geneID=gene)
    if(nrow(exp)>1){
      rownames(exp) <- byTF$TF
    }
    P <- scatter(meth=meth, exp=exp,category=category,
                 xlab="Avg DNA methyation", ylab=sprintf("TF expression"), 
                 title="TF vs avg DNA methylation",
                 ...)
    if(save) ggsave(filename = sprintf("%s/%s.byTF.pdf",dir.out,paste(byTF$TF,collapse = "_")),
                    plot = P,useDingbats=FALSE, width=3*(length(byTF$TF)%%5), 
                    height = 3*ceiling(length(byTF$TF)/5))
  }
  return(P)
}


#'scatter
#'@importFrom reshape melt.data.frame
#'@importFrom ggplot2 geom_point facet_wrap scale_x_continuous theme_bw theme element_blank labs scale_colour_manual geom_smooth element_text
#@param meth A vector of number.
#@param exp A vector of number or matrix with sample in column and gene in rows.
#@param category A vector of sample labels.
#@param xlab A character specify the title of x axis.
#@param ylab A character specify the title of y axis.
#@param title A character specify the figure title.
#@param color.value A vector specify the color of each category, such as 
#color.value=c("Experiment"="red","Control"="darkgreen")
#@param lm_line A logic. If it is TRUE, regression line will be added to the graph.
#@return A ggplot figure object
scatter <- function(meth, exp, category=NULL, xlab=NULL, ylab=NULL,title=NULL,
                    color.value=NULL,lm_line=FALSE){
  if(is.null(category)) category <- rep(1,length(meth))
  if(!is.vector(exp)){
    exp <- as.data.frame(t(exp))
    GeneID <- colnames(exp)
    exp$meth <- meth
    exp$category <- category
    df <- melt.data.frame(exp, measure.vars = GeneID)
    P <- ggplot(df, aes(x= meth, y=value, color=factor(category)))+
      geom_point(size=0.9)+
      facet_wrap(facets = ~ variable, ncol = 5)+
      scale_x_continuous(limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1))+
      theme_bw()+
      theme(panel.grid.major = element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
      labs(x=xlab,y=ylab,title=title)
    if(!is.null(color.value)) P <- P+scale_colour_manual(values = color.value)
    if(lm_line) P <- P+geom_smooth(method = "lm", se=FALSE, color="black", 
                                   formula = y ~ x,data=df)
  }else{
    df <- data.frame(meth=meth,exp=exp,category=category)
    if(length(unique(df$category))==1){
      P <- ggplot(df, aes(x= meth, y=exp))
    }else{
      P <- ggplot(df, aes(x= meth, y=exp, color=factor(category)))
    }
      P <- P + geom_point()+
      scale_x_continuous(limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
      labs(x=xlab,y=ylab,title=title)
    if(lm_line){
#       P <- P+ geom_text(aes(x =0.8 , y = max(exp)-0.5, label = lm_eqn(df)),
      #parse = TRUE,colour = "black")+
        P <- P+ geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x,data=df)
    }
  }
  return(P)
}

