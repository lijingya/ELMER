#' Read a bed file.
#' @param x A path of bed file (characters)
#' @param strand A boolean to specific strands. If true, strand column will be filled as input. 
#' If false, strand column will be filled "*""
#' @param skip A number to specify how many lines should be removed from bed file.
#' @param cols Specify the column to read from bed file. 
#' @param seqLength Specify custmer seqLength parameter in GRange function
#' @return GRange object containning bed file information.
#' @export



ReadBed <- function(x,strand=FALSE,skip=0,cols=NULL,seqLength=NULL){
  if(is.null(cols)){
    x <- read.table(x,stringsAsFactors=FALSE,sep = "\t",skip=skip)
  }else{
    cols <- paste(cols,collapse=",")
    cmd <- sprintf("cut -f%s %s",cols,x)
    x <- read.table(pipe(cmd),stringsAsFactors=FALSE,sep = "\t",skip=skip)
  }
  
  x[,1] <- sub("chrx","chrX",x[,1])
  x[,1] <- sub("chry","chrY",x[,1])
  x[,1] <- sub("chrm","chrM",x[,1])
  if(dim(x)[2] >5){
    x[,6] <- sub("\\.","*",x[,6])
  }
  x[,1] <- sub(" ","",x[,1])
  if(strand){
    Bed<-GRanges( x[,1], IRanges(x[,2]+1, x[,3]) ,strand=x[,6])
    if(ncol(x)>6){
      values(Bed) <- data.frame(name=x[,4],score = x[,5],
                                x[,setdiff(seq_len(ncol(x)),c(1:6))])
    }  
  }else{
    Bed<-GRanges( x[,1], IRanges(x[,2], x[,3]) )
    if(ncol(x)==4){
      values(Bed) <- data.frame(name=x[,4],stringsAsFactors = F)
    }else if(ncol(x)==5){
      values(Bed) <- data.frame(name=x[,4],score = x[,5],stringsAsFactors = F)
    }else if(ncol(x)>=6){
      values(Bed) <- data.frame(name=x[,4],score = x[,5],
                                x[,setdiff(seq_len(ncol(x)),c(1:6))],
                                stringsAsFactors = F)
    }  
  }
  if(!is.null(seqLength)){
    seqlengths(Bed) <- sequenceLen[seqlevels(Bed),2]
  }
  
  return(Bed)
}

#' Read a GFF file.
#' @param x A path of GFF file (characters)
#' @param strand A boolean to specific strands. If true, strand column will be 
#' filled as input. If false, strand column will be filled "*""
#' @param skip A number to specify how many lines should be removed from bed file.
#' @return GRange object containning GFF file information.
#' @export
ReadGFF <- function(x,strand=FALSE,skip=0  ){
  x <- read.table(x,stringsAsFactors=FALSE,sep = "\t",skip=skip)
  x[,1] <- sub("chrx","chrX",x[,1])
  x[,1] <- sub("chry","chrY",x[,1])
  x[,1] <- sub("chrm","chrM",x[,1])
  x[,7] <- sub("\\.","*",x[,7])
  x[,1] <- sub(" ","",x[,1])
  if(strand){
    GFF<-GRanges( x[,1], IRanges(x[,4], x[,5]) ,strand=x[,7])
    values(GFF) <- data.frame(Source=x[,2],feature=x[,3],score = x[,6],
                              frame = x[,8],name = x[,9],
                              x[,setdiff(seq_len(ncol(x)),c(1:9))],
                              stringsAsFactors = F)
  }else{
    GFF<-GRanges( x[,1], IRanges(x[,4], x[,5]) )
    values(GFF) <- data.frame(Source=x[,2],feature=x[,3],score = x[,6],
                              frame = x[,8],name = x[,9],
                              x[,setdiff(seq_len(ncol(x)),c(1:9))],
                              stringsAsFactors = F)
  }
  seqlengths(GFF) <- sequenceLen[seqlevels(GFF),2]
  return(GFF)
}

sequenceLen <- 
  data.frame(chr=c("chr1","chr1_gl000191_random","chr1_gl000192_random","chr2",
                   "chr3","chr4","chr4_gl000193_random","chr4_gl000194_random","chr5",
                   "chr6","chr7","chr7_gl000195_random", "chr8","chr8_gl000196_random",
                   "chr8_gl000197_random","chr9","chr9_gl000198_random",
                   "chr9_gl000199_random","chr9_gl000200_random","chr9_gl000201_random",
                   "chr10","chr11","chr11_gl000202_random","chr12","chr13","chr14",
                   "chr15","chr16","chr17","chr17_gl000203_random","chr17_gl000204_random",
                   "chr17_gl000205_random", "chr17_gl000206_random","chr18",
                   "chr18_gl000207_random","chr19","chr19_gl000208_random","chr19_gl000209_random",
                   "chr20", "chr21", "chr21_gl000210_random", "chr22","chrX",
                   "chrY","chrM", "chrUn_gl000211","chrUn_gl000212","chrUn_gl000213",
                   "chrUn_gl000214","chrUn_gl000215","chrUn_gl000216","chrUn_gl000217",
                   "chrUn_gl000218","chrUn_gl000219","chrUn_gl000220",
                   "chrUn_gl000221","chrUn_gl000222","chrUn_gl000223","chrUn_gl000224",
                   "chrUn_gl000225", "chrUn_gl000226", "chrUn_gl000227",
                   "chrUn_gl000228","chrUn_gl000229","chrUn_gl000230","chrUn_gl000231",
                   "chrUn_gl000232","chrUn_gl000233","chrUn_gl000234",
                   "chrUn_gl000235","chrUn_gl000236","chrUn_gl000237","chrUn_gl000238",
                   "chrUn_gl000239","chrUn_gl000240",
                   "chrUn_gl000241", "chrUn_gl000242", "chrUn_gl000243","chrUn_gl000244",
                   "chrUn_gl000245","chrUn_gl000246","chrUn_gl000247",
                   "chrUn_gl000248","chrUn_gl000249"  ),
             Length=c(249250621,106433,547496,243199373,198022430,191154276, 
                      189789,191469,180915260,171115067,159138663,182896,146364022,
                      38914,37175,141213431,90085,169874, 187035, 36148,135534747,
                      135006516,40103,133851895,115169878,107349540,102531392,
                      90354753,81195210,37498,81310, 174588,41001,78077248, 4262, 
                      59128983,92689,159169,63025520,48129895,27682,51304566,
                      155270560,59373566,16569,166566,186858,164239,137718,172545,
                      172294,172149,161147,179198,161802,155397,186861, 180455,
                      179693,211173,15008, 128374,129120,19913,43691,27386,40652,
                      45941, 40531,34474,41934, 45867, 39939,33824,41933,42152,
                      43523,43341,39929,36651,38154,36422,39786,38502
             ))
                            
                         
rownames(sequenceLen) <- sequenceLen$chr

#' Write a bed file from GRange object.
#' @param x GRange object
#' @param save if save is false, function will return a bed format data.frame. 
#' If save is true, fn parameter need to be specific and it output bed file in the path you specified in fn.
#' @param fn A name of bed file you want to output.
#' @return A data.frame bed object or save output bed file.
#' @export
WriteBed <-function(x,save=T,fn=NULL){
  x <- as.data.frame(x,row.names=NULL)
  x$element=NULL
  if(ncol(x) >5){
    if(ncol(x)>=7){
      out <- x[,c(1,2,3,6,4,5,7:ncol(x))]
    }else{
      out <- x[,c(1,2,3,6,4,5)]
    }
  }else{
    Names <- paste0(x[,1],":",x[,2],"-",x[,3])
    out <- cbind(x[,c(1,2,3)],Names,x[,4],x[,5])
  }
  rownames(out) <- NULL
  if(save){
    if(is.null(fn)) fn <- deparse(substitute(x))
    write.table(out,file=fn,row.names=F,col.names=F,quote=F,sep="\t")
  }else{
    return(out)
  }
}

#' Generate random loci of genome.
#' @param SampleSize A number of random loci you want to generate.
#' @param exclusion The chromosome you want to exclude such as chrX chrY.
#' @param regionWidth The width of each random loci.
#' @return GRange object.
RandomLoci <- function(SampleSize=NULL,exclusion=NULL,regionWidth=0){
  chr <- paste0("chr",c(1:22,"X","Y"))
  if(!is.null(exclusion)) chr <- chr[!chr %in% exclusion]
  LengthSum <- sum(sequenceLen[chr,2])
  positions <- sample(seq_len(LengthSum),SampleSize)
  out <- c()
  for(i in seq_len(SampleSize)){
    count <- 1
    start <- positions[i]
    while((start- sequenceLen[chr[count],2]) > 0){
      start <- start- sequenceLen[chr[count],2]
      count <- count+1
    }
    chrsample <- chr[count]
    end <- start +regionWidth
    if(end > sequenceLen[chrsample,2]){
      end <- sequenceLen[chrsample,2]
      start <- end - regionWidth
    }
    out <- rbind(out,c(chrsample,start,end))
  }
  out <- GRanges( out[,1], IRanges(as.numeric(out[,2]), as.numeric(out[,3])))
  seqlengths(out) <- sequenceLen[seqlevels(out),2]
  return(out)
}


