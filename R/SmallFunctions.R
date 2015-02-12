# #extract probe bed file: need to write code for this later. out put a bed format file. Next step can direct load bed file as Granges.
# 
# #2. Form matrix for feature and probe overlap
# # 1 means overlap, 0 means no overlap
# # target_region and bioFeature_fl need to be bed file format if it is character. If it is already loaded file, it need to be Grange format
# 
# 
# 
# 
# feature_overlap <- function (target_region, bioFeature_fl,target.size=c(0,0),bio.size=c(400,400),col.names=NULL, row.names=NULL){  # bioFeature_fl : list of pathes of a banch of files that will otherlap with target region. 
#   require(GenomicRanges)
#   n_biofeature <- length (bioFeature_fl)
#   
#   #####load target region information
#   if (mode(target_region)== "character") RegionInfo <- ReadBed(target_region)
#   start(RegionInfo) <- start(RegionInfo)-target.size[1]
#   end(RegionInfo) <- end(RegionInfo)+target.size[2]
#   #####load biofeature information and form biofeature matrix  
#   
#   biofeature_matrix = c()
#   colName <- c()
#   for (i in 1:n_biofeature){
#     if (any(mode(bioFeature_fl)== "character")){
#       if(file.info(bioFeature_fl[i])$size==0){
#         print (paste0("Error: no lines available in input",bioFeature_fl[i]))                   #read.table can't read 0 size input.
#       }else{
#         file_regionINFO = ReadBed(bioFeature_fl[i])
#         start(file_regionINFO) <- start(file_regionINFO)-bio.size[1]
#         end(file_regionINFO) <- end(file_regionINFO)+bio.size[2]
#         overlap_vector = RegionInfo %over% file_regionINFO
#         biofeature_num = sum(file_regionINFO %over% RegionInfo)
#         biofeature_matrix = cbind(biofeature_matrix,overlap_vector)
#         print(paste0("reading and overlap ", bioFeature_fl[i]))
#         bn <- basename(bioFeature_fl[i])
#         bn<-sub(".bed","",bn)
#         colName <- c(colName,bn)
#       }
#     }else{
#       start(file_regionINFO) <- start(file_regionINFO)-bio.size[1]
#       end(file_regionINFO) <- end(file_regionINFO)+bio.size[2]
#       overlap_vector = RegionInfo %over% file_regionINFO
#       biofeature_num = sum(file_regionINFO %over% RegionInfo)
#       biofeature_matrix = cbind(biofeature_matrix,overlap_vector)
#       print(paste0("reading and overlap ", bioFeature_fl[i]))
#       colName <- c(colName,paste0(bioFeature_fl[i],"_",biofeature_num))
#     }
#     
#   }
#   
#   biofeature_matrix[biofeature_matrix==FALSE]<- 0
#   biofeature_matrix[biofeature_matrix==TRUE]<- 1
#   
#   ####merge the target region with the biofeature matrix
#   
#   colnames(biofeature_matrix) <- colName
#   rownames(biofeature_matrix) <- RegionInfo$name
#   final_file = as.data.frame (final_file)                                                  # it is better to use dataframe, it has dim. the list don't have dim
#   return (final_file)
# }

# 
# #------------------------------MergeBioFeature function----------------------
# #main is the main manifest
# #custom is list of data which will be merged into main manifest
# #by is specify the column to use to merge the files, default is the row.names.
# MergeBioFeature <- function(main, custom=list(),by.main="row.names",by.custom=c()){
#   if(!is.data.frame(main)){
#     main = as.data.frame (main)
#   }
#   
#   for (i in 1:length(custom)){    
#     oldn<-nrow(main)
#     main<-merge(main,custom[[i]],by.x=by.main,by.y=by.custom[i],sort=TRUE)
#     if(by.main == "row.names" && by.custom[i]=="row.names"&& i < length(custom)){
#       rownames(main)<-main$Row.names # I don't know why merge doesn't do this by default
#       main$Row.names<-NULL
#     }
#     print(sprintf("Merge of main file (%d rows) and extra cols (%d rows) yielded %d rows\n",oldn,nrow(custom[[i]]),nrow(main)))
#   }
#   
#   return(main)  
# }


#-------------------------------repliAnalysis function(from Ben)---------------------

repliAnalysis<-function() {
  
  repliFns <- dir(dir1,pattern="Rep1.bed$",full.names=TRUE)
  print(paste0("Found ",length(repliFns)," repli files"))
  
  round <- 1
  
  colName <- c("id")                                                # record the column names
  out <- c()                                                        # output
  
  for (repliFn in repliFns)
  {
    
    bn <- basename(repliFn)
    bn<-sub("HM450-wgEncodeUwRepliSeq","",bn)
    bn<-sub(".bed","",bn)
    bn<-sub("WaveSignalRep1","",bn)  
    colName <- c(colName,bn) 
    
    repli<-read.table(repliFn,sep="\t")
    if(round==1){
      id <- repli[,4]
      out <- as.character(id)
    }else{
      id.new <- repli[,4]
      stopifnot(length(id)==length(id.new) && all(as.character(id)==as.character(id.new)))
    }
    print(paste0("RepliAnalysis saw file: ",repliFn))
    out<-cbind(out, repli[,5])
    print(paste0("output has ",ncol(out)-1," repli cell type: ",bn))
    round <- round + 1
  }
  
  colnames(out)<-colName
  return(out)
}

#-----------------------combineMethysets function----------------------------
#ProbesToUse is probes names vector to indicate which should to use
# CancerTypeToDo indices are from the array of tumor types
#SampleToUse : if it is null , it will output all sampletypes, if it is TN, it means Tumor and tissue normal sample.
combineMethysets <- function(CancerTypeToDo=NULL,Probes=NULL, SampleToUse=NULL,panCan=panCan)
{
  load(file11)
  if(is.null(Probes)) {
    Probes <- read.table(pipe(sprintf("cut -f4 %s",file2.2)),colClasses = "character")[,1]
  }
  if(!is.null(CancerTypeToDo)){
    if(all(CancerTypeToDo %in% c("BLCA","BRCA","COAD","GBM","HNSC","KIRC","LAML","LUAD","LUSC","READ","UCEC"))){
      CancerTypeToDo <- c(1:11)[c("BLCA","BRCA","COAD","GBM","HNSC","KIRC","LAML","LUAD","LUSC","READ","UCEC") %in% CancerTypeToDo]
    }
    subAnnot <- subAnnot[subAnnot$CT %in% CancerTypeToDo,]
  } 
  if (SampleToUse=="TN"){
    subAnnot <- subAnnot[subAnnot$TN.cat %in% "Tumor" |subAnnot$TN.cat %in% "Normal",]
  }else if(SampleToUse=="T"){
    subAnnot <- subAnnot[subAnnot$TN.cat %in% "Tumor",]
  }else if(SampleToUse=="N"){
    subAnnot <- subAnnot[subAnnot$TN.cat %in% "Normal",]
  }
  betas <- panCan[Probes,subAnnot$sampleNames]
  final <-list(betas=betas,Info=subAnnot)
  rm(subAnnot)
  return(final)
}


##------------mutation Info----------------------------------------------------------
#SampleSet is the sample barcode with proper order
GetMutation<-function(mutation=file6,SampleSet=NULL,MutNames=NULL) 
{
  require(snow)
  newenv <- new.env()
  load(mutation, envir=newenv)
  Mut <- get(ls(newenv)[1],envir=newenv) # The data is in the one and only variable
  Mut$Tumor_Sample_Barcode <- standardizeTcgaIds(Mut$Tumor_Sample_Barcode)
  Mut$Matched_Norm_Sample_Barcode <- standardizeTcgaIds(Mut$Matched_Norm_Sample_Barcode)
  Mut <- Mut[,c(1,16,9,10)]
  SampelOverlap <- Mut$Tumor_Sample_Barcode %in% SampleSet
  print (paste0("The number of the smaple that has mutation data is ",sum(SampleSet %in% Mut$Tumor_Sample_Barcode)))
  if(is.null(MutNames)){
    MutNames <- names(sort(table(Mut[ Mut$Tumor_Sample_Barcode %in% SampleSet,1])[table(Mut[ Mut$Tumor_Sample_Barcode %in% SampleSet,1])>0]))
  }else{
    MutNames <- MutNames
  }
  
  # Form the matrix of the mutation
  pthreads <- detectCores()/2
  cl <- makeCluster(pthreads,type="SOCK")
  out <- parLapplyLB(cl,MutNames,Mutvector,SampleSet,Mut)
  stopCluster(cl)
  out <- do.call(cbind,out)
  colnames(out) <- paste0(MutNames,"mut")
  rownames(out) <- SampleSet
  return(out);
}

Mutvector <- function(MutName,SampleSet,Mut){
 tmp <- Mut[as.character(Mut$Hugo_Symbol) %in% MutName,]
 out<- as.character(tmp$Variant_Classification[match(SampleSet,tmp$Tumor_Sample_Barcode)])
 return(out)
}


##-----------RNA-seq-----------------------------------------------------------------

#Gene make rowname separat -------------------------------------
GeneIDName <- function(x){
  tmp<-strsplit(rownames(x),"\\|")
  GeneID<-unlist(lapply(tmp,function(x) x[2]))
  GeneID <- paste0("ID",GeneID)
  row.names(x) <- GeneID
  return(x)
}


##-------------Get gene symbol through ID-------------------------------
.GENEID2Symbol <- function(x,revert=FALSE){
  load("/export/uec-gs1/laird/users/lijingya/data/methylation/TCGA/array/450K/INFO/UCSC.hg19.knownGene_TSS_gene.rda")
  TxDbTSS$GENEID <- paste0("ID",TxDbTSS$GENEID)
  Pairs <- unique(data.frame(GENEID=TxDbTSS$GENEID,Symbol= TxDbTSS$SYMBOL))
  Pairs[,1] <- as.character(Pairs[,1])
  Pairs[,2] <- as.character(Pairs[,2])
  if(revert){
    out <- Pairs[match(x,Pairs$Symbol), "GENEID"]
    out <- x[is.na(out)]
  }else{
    out <- Pairs[match(x,Pairs$GENEID), "Symbol"]
    out[is.na(out)] <- x[is.na(out)]
  }
  return(out)
}

##-------------make TxDBTSS file-------------------------------------
MakeTxDbTSS <- function(){
  library(Homo.sapiens)
  keytypes(Homo.sapiens)
  txs <- transcriptsBy(Homo.sapiens, 'gene', col=c('GENEID','SYMBOL'))
  txs <- unlist(txs)
  TxDBTSS <- promoters(txs,0,0)
  save(TxDBTSS,file="/export/uec-gs1/laird/users/lijingya/data/methylation/TCGA/array/450K/INFO/UCSC.hg19.knownGene_TSS_gene.rda")
  #in pipeline it should be storage to somewhere.
}


# ##------------z score calculation------------------------------------
# # reduce memory occupy I should remove Exps and Meths.
# #By: "gene" or "probe" which must be one.
# .Stat.zscore <- function(Probe,Gene,K,Top=NULL){
#   if(!(length(Gene)==1 & length(Probe)==1)) {stop("Number of Gene ID or Probe should be 1")}
#   Exp <- as.matrix(Exps[Gene,])[1,]
#   Meth <- Meths[Probe,]
#   Meth_B <- Binary(Meth,Break=K)
#   zscore <- c()
#   if(sum(Exp)==0){
#     zscore <- NA
#   }else{
#     df <- data.frame(Exp=Exp,Meth=Meth,Meth_B=Meth_B)
#     rownames(df) <- names(Meth)
#     unmethyPercent <- sum(df$Meth_B==0,na.rm=T)/dim(df)[1]
#     methyPercent <- sum(df$Meth_B==1,na.rm=T)/dim(df)[1]
#     if(unmethyPercent < 0.05 | methyPercent < 0.05){
#       zscore <- NA
#     }else{
#       unmethy <- rownames(df[order(df$Meth),])[1:round(nrow(df)*Top)] 
#       methy <- rownames(df[order(df$Meth,decreasing=T),])[1:round(nrow(df)*Top)] 
#       zscore <- (mean(df[unmethy,"Exp"],na.rm=T)-mean(df[methy,"Exp"],na.rm=T))/sd(df[methy,"Exp"],na.rm=T)
#     }
#     out <- c(Gene, Probe,zscore)
#   } 
#   return(out)
# }
# 
# .Stat.zscore2 <- function(Probe,Gene,K,Top=NULL){
#   if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
#   Exp <- as.matrix(Exps[Gene,])
#   Meth <- Meths[Probe,]
#   Meth_B <- Binary(Meth,Break=K)
#   zscore <- c()
#   if(sum(Exp)==0){
#     zscore <- NA
#   }else{
#     unmethyPercent <- sum(Meth_B==0,na.rm=T)/length(Meth_B)
#     methyPercent <- sum(Meth_B==1,na.rm=T)/length(Meth_B)
#     if(unmethyPercent < 0.05 | methyPercent < 0.05){
#       zscore <- NA
#     }else{
#       unmethy <- names(Meth[order(Meth)])[1:round(length(Meth)*Top)] 
#       methy <- names(Meth[order(Meth,decreasing=T)])[1:round(length(Meth)*Top)] 
#       zscore <- apply(Exp,1,function(x) {(mean(x[unmethy],na.rm=T)-mean(x[methy],na.rm=T))/sd(x[methy],na.rm=T)})
#     }
#     out <- cbind(Probe=rep(Probe,length(Gene)),Gene,zscore)
#   } 
#   return(out)
# }
# 
# .Stat.zscore3 <- function(Probe,NearGenes,K,Top=NULL){
#   if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
#   Gene <- NearGenes[[Probe]][,2]
#   Exp <- as.matrix(Exps[Gene,])
#   Meth <- Meths[Probe,]
#   Meth_B <- Binary(Meth,Break=K)
#   zscore <- c()
#   if(sum(Exp)==0){
#     zscore <- NA
#   }else{
#     unmethyPercent <- sum(Meth_B==0,na.rm=T)/length(Meth_B)
#     methyPercent <- sum(Meth_B==1,na.rm=T)/length(Meth_B)
#     if(unmethyPercent < 0.05 | methyPercent < 0.05){
#       zscore <- NA
#     }else{
#       unmethy <- names(Meth[order(Meth)])[1:round(length(Meth)*Top)] 
#       methy <- names(Meth[order(Meth,decreasing=T)])[1:round(length(Meth)*Top)] 
#       zscore <- apply(Exp,1,function(x) {(mean(x[unmethy],na.rm=T)-mean(x[methy],na.rm=T))/sd(x[methy],na.rm=T)})
#     }
#     out <- cbind(Probe=rep(Probe,length(Gene)),Gene,zscore)
#   } 
#   return(out)
# }
# 
# 
# 
# 
# .Stat.corr2 <- function(Probe,Gene,K,Top=NULL,method="pearson"){
#   if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
#   Exp <- as.matrix(Exps[Gene,])
#   Meth <- Meths[Probe,]
#   Meth_B <- Binary(Meth,Break=K)
#   corr <- c()
#   if(sum(Exp)==0){
#     corr <- NA
#   }else{
#     unmethyPercent <- sum(Meth_B==0,na.rm=T)/length(Meth_B)
#     methyPercent <- sum(Meth_B==1,na.rm=T)/length(Meth_B)
#     if(unmethyPercent < 0.05 | methyPercent < 0.05){
#       corr <- NA
#     }else{
#       corr <- apply(Exp,1,function(x,Meth) {cor(x,Meth,use="complete.obs",method = method)},Meth=Meth)
#     }
#     out <- cbind(Probe=rep(Probe,length(Gene)),Gene,corr)
#   } 
#   return(out)
# }
# 
# .Stat.corr3 <- function(Probe,NearGenes,K,Top=NULL,method="pearson"){
#   if(! length(Probe)==1) {stop("Number of  Probe should be 1")}
#   Gene <- NearGenes[[Probe]][,2]
#   Exp <- as.matrix(Exps[Gene,])
#   Meth <- Meths[Probe,]
#   Meth_B <- Binary(Meth,Break=K)
#   corr <- c()
#   if(sum(Exp)==0){
#     corr <- NA
#   }else{
#     unmethyPercent <- sum(Meth_B==0,na.rm=T)/length(Meth_B)
#     methyPercent <- sum(Meth_B==1,na.rm=T)/length(Meth_B)
#     if(unmethyPercent < 0.05 | methyPercent < 0.05){
#       corr <- NA
#     }else{
#       corr <- apply(Exp,1,function(x,Meth) {cor(x,Meth,use="complete.obs",method = method)},Meth=Meth)
#     }
#     out <- cbind(Probe=rep(Probe,length(Gene)),Gene,corr)
#   } 
#   return(out)
# }


#---probes name index------------------------------------
# calculate Pvalue
Get.Pvalue <- function(zscore.Matrix,permu){
  .Pvalue <- function(x,permu,target.Matrix){
    zscore <- target.Matrix[x,3]
    Gene <- target.Matrix[x,"Gene"]
    if(is.na(zscore)){
      out <- NA
      print("NA")
    }else{
      out <- sum(permu[Gene,]  > zscore | permu[Gene,] == zscore,na.rm=T)/sum(!is.na(permu[Gene,]))
    }
#  else if(zscore <= 0){
#    out <- sum(permu[Gene,] < zscore | permu[Gene,] == zscore, na.rm=T)/length(permu[Gene,])   
    return(out)
  }
  Pvalue <- apply(matrix(1:nrow(zscore.Matrix),ncol=1),1,.Pvalue,target.Matrix=zscore.Matrix,permu=permu)
  Pvalue <- unlist(Pvalue)
  Output <- cbind(zscore.Matrix,Pvalue)
  return(Output)
}



Get.Pvalue.c <- function(zscore.Matrix,permu){
  .Pvalue <- function(x,permu,target.Matrix){
    zscore <- target.Matrix[x,3]
    Gene <- target.Matrix[x,"Gene"]
    if(is.na(zscore)){
      out <- NA
      #print("NA")
    }else{
      out <- sum(permu[Gene,]  < zscore | permu[Gene,] == zscore,na.rm=T)/sum(!is.na(permu[Gene,]))
    }
    #  else if(zscore <= 0){
    #    out <- sum(permu[Gene,] < zscore | permu[Gene,] == zscore, na.rm=T)/length(permu[Gene,])   
    return(out)
  }
  Pvalue <- apply(matrix(1:nrow(zscore.Matrix),ncol=1),1,.Pvalue,target.Matrix=zscore.Matrix,permu=permu)
  Pvalue <- unlist(Pvalue)
  Output <- cbind(zscore.Matrix,Pvalue)
  return(Output)
}
z.test = function(a, mu, var){
  zeta = (mean(a) - mu) / (sqrt(var / length(a)))
  return(zeta)
}

CNV_forgene <- function(wd,CT,Type){
  files <- dir(wd)
  CT <- toupper(CT)
  if(Type %in% "thresholded"){
    File <- paste0(wd,"/",CT,".gistic.all_thresholded.by_genes")
  }else{
    File <- paste0(wd,"/",CT,".gistic.all_data_by_genes")
  }
  CNV_Gene <- read.table(File,header= T,sep="\t",stringsAsFactors=F)
  colnames(CNV_Gene) <- gsub("\\.", "-", colnames(CNV_Gene))
  rownames(CNV_Gene) <- CNV_Gene[,1]  
  CNV_Gene <- CNV_Gene[,-c(1:3)]
  return(CNV_Gene)
}

##CordinateToRange
CordinateToRange <- function(Cordinate,core=detectCores()/2){
  options('mc.cores'=core)
  chr <- unlist(mclapply(Cordinate,function(x){strsplit(x,"\\:")[[1]][1]}))
  position <- sub("chr.*:","",Cordinate)
  startAend <- matrix(unlist(mclapply(position,function(x){strsplit(x,"\\-")})),ncol=2,byrow=T)
  class(startAend) <- "numeric"
  Bed<-GRanges(chr, IRanges(startAend[,1], startAend[,2]) )
  return(Bed)
}


#---------------------------------------lsos Function---------
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(print(object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

MakeHeatmap <- function(sub.probes,METH,CT,K=0.3,title=NULL){
  # set colors --------------------------------------------------------------
  jet.colors <-
    colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  redGreen <- colorRampPalette(c("green","black","red"))
  library(RColorBrewer)
  GreyCol <- colorRampPalette(brewer.pal(9,"Greys"))
  # load data ---------------------------------------------------------------
  ## load REMC barcode 
  #CT specific barcode
  REMC_Barcode<- c("GI","GI","BRST","BRN","THYM","KID","BLD","LNG","LNG","GI","OVRY")
  
  REMC_ID <- read.csv("/export/uec-gs1/laird/users/lijingya/data/REMC/FinalFiles_Mar2014/barcode.csv",stringsAsFactors=F)
  relativeID <- unique(unlist(mclapply(REMC_Barcode[CT],function(x){ REMC_ID$Epigenome.Mnemonic[grepl( x, REMC_ID$Epigenome.Mnemonic)]})))
  ESC <- REMC_ID[REMC_ID$GROUP %in% "ESC" | REMC_ID$Epigenome.Mnemonic %in% relativeID,]
  # load probe info
  tmp <- as.matrix(read.csv("/export/uec-gs1/laird/users/lijingya/analysis/methylation/TCGA/450Kmethyl/Results/2014-04-23_RMEC_Phantom5/Probes_Phantom_REMC_TSS_CGI_info_hg19.csv",nrow=1,header=F))[1,]
  
  Cols <- paste(c(1,2,3,909:920,which(tmp %in% c("E119",ESC$NEW.EID))),collapse=',')
  Probes <- read.csv(
    pipe(sprintf('cut -f %s -d "," /export/uec-gs1/laird/users/lijingya/analysis/methylation/TCGA/450Kmethyl/Results/2014-04-23_RMEC_Phantom5/Probes_Phantom_REMC_TSS_CGI_info_hg19.csv',Cols)),
    stringsAsFactors=F)
  Probes <- Probes[rowSums(Probes[,c("EnhG1","EnhG2","EnhA1","EnhA2")])>0,]
  ESC <- ESC[ESC$NEW.EID %in% colnames(Probes),]
  Probes <- Probes[match(sub.probes,Probes$Probes),]
  
  ## read buffy coat data
  Buffy <- read.table("/home/rcf-40/lijingya/lijingya/data/methylation/TCGA/array/450K/INFO/1074_Buffy\ Coat\ corrected\ beta\ values.txt",skip=5,stringsAsFactors=F,header=F)
  Buffy <- as.matrix(Buffy)
  rownames(Buffy)<- Buffy[,1]
  Buffy <- Buffy[,-1]
  class(Buffy) <- "numeric"
  Buffy <- Buffy[sub.probes,]
  
  ## fimo motif
  ##fimo motif search
  load("/export/uec-gs1/laird/users/lijingya/analysis/methylation/TCGA/450Kmethyl/Results/2014-05-06-BCRAexample/exploreClustering/Probesall.TF.matrix.rda")
  
  ## REMC WGBS
  load("/export/uec-gs1/laird/users/lijingya/data/methylation/TCGA/array/450K/INFO/Probes_REMC_WGBS.rda")
  
  #load lower.bound
  load(sprintf("../CT_%d_low.assign_SigHypo.rda",CT))
  
  sub.METH <- METH$betas[sub.probes,]
  sub.Info <- METH$Info
  sub.METH.normal <- sub.METH[,sub.Info$TN.cat %in% "Normal"]
  sub.METH.tumor <- sub.METH[,sub.Info$TN.cat %in% "Tumor"]
  sub.Info <- sub.Info[sub.Info$TN.cat %in% "Tumor",]   
  sub.normal <- Binary(sub.METH.normal,Break=K)
  sub.tumor <- Binary(sub.METH.tumor,Break=K)
  system.time(cluster_T <- cluster.main(sub.tumor,Colv=T,Rowv=T,hclustMethod="ward",distMethod="binary"))
  cluster_N <- cluster.main(sub.normal,Colv=T,Rowv=F,hclustMethod="ward")
  cluster_T$x <- sub.METH.tumor
  cluster_N$x <- sub.METH.normal
  cluster_N$rowInd <- cluster_T$rowInd
  ##
  sub.CGI <- Probes[match(sub.probes,Probes$Probes),"CGI"]
  sub.buffy <- Buffy[rownames(sub.METH),]
  sub.segmentation <- Probes[match(sub.probes,Probes$Probes),ESC$NEW.EID]
  tmp <- ChrToColor(sub.segmentation)
  Leve <- tmp$levels
  sub.segmentation <- tmp$x
  sub.REMC.WGBS <- Probes.REMC.WGBS[rownames(sub.METH),]
  colnames(sub.REMC.WGBS) <- REMC_ID[match(colnames(sub.REMC.WGBS),REMC_ID$NEW.EID),3]
  sub.TF.matrix <- Probes.TF[sub.probes,]
  sub.enrich <- unlist(mclapply(1:ncol(sub.TF.matrix),function(x) {binom.test(colSums(sub.TF.matrix)[x], nrow(sub.TF.matrix), p = colMeans(Probes.TF)[x],alternative =  "greater")$p.value},mc.cores=6))
  sub.TF.matrix <- sub.TF.matrix[,names(sort(-log10(sub.enrich),decreasing=T))[1:20]]
  sub.probes <- rownames(sub.METH.tumor)[unlist(apply(sub.METH.tumor,1,function(x){ all(!is.na(x))}))]
  Mean.TF.sample <- t(sub.METH.tumor[sub.probes,]) %*% sub.TF.matrix[sub.probes,] / matrix(rep(colSums(sub.TF.matrix[sub.probes,]),ncol(sub.METH.tumor)),nrow=ncol(sub.METH.tumor),byrow=T)
  TF.cluster <- cluster.main(sub.TF.matrix,Colv=T,Rowv=F,hclustMethod="ward")
  sub.TF.matrix <- sub.TF.matrix[,TF.cluster$colInd]
  ImpInfo <- c("VHLmeth" ,"CDKN2Ameth" ,"BRCAmeth", "MLH1meth","VHL", "BRCA1","BRCA2", "CDKN2A","TP53","ARID1A","MLL2","SETD2","MLL3","CTCF", "KDM6A","MLH1","ATRX","CREBBP","EP300","DNMT3A",
               "TET2","TET1","IDH1","IDH2","KDM3A","KDM5C","MLL4","MLL","ASXL1","EZH2","MECOM","DAXX","IDH1.R132","IDH2.R172","IDH2.R140","BRAFV600E","KRAS")
  sub.ImpInfo <- sub.Info[,ImpInfo]
  sub.ImpInfo <- do.call(cbind,mclapply(colnames(sub.ImpInfo),function(x,sub.ImpInfo){ as.numeric(sub.ImpInfo[,x])},sub.ImpInfo=sub.ImpInfo,mc.cores=6))
  sub.ImpInfo[is.na(sub.ImpInfo)] <- 0
  colnames(sub.ImpInfo) <- ImpInfo
  if(paste(CT,collapse="_") %in% "2"){
    sub.pam <- as.numeric(sub.Info$pam50)
    sub.pam[is.na(sub.pam)] <- -1
  }
  if(length(CT) > 1){
    sub.pam <- sub.Info$CT
  }
  png(paste0("CT",paste(CT,collapse="_"),"_",title,"_probes",dim(sub.tumor)[1],".png"),  width = 1600, height = 1500)
  if(paste(CT,collapse="_") %in%"2" | length(CT) > 1){
    layout(matrix(data=c(0,16,17,17,17,17,0,
                         10,9,0,0,0,0,0,
                         8,7,0,0,0,0,0,
                         6,5,14,14,14,14,0,
                         4,3,0,0,0,0,0,
                         2,1,11,12,13,15,18,
                         0,19,0,0,0,0,0,
                         0,20,0,0,0,0,0), nrow=8,ncol=7,byrow=TRUE),widths=c(5*(ncol(sub.normal)/ncol(sub.tumor)),5,0.2,0.4,2,2.8,2.8),heights=c(0.2,0.2,0.2,0.2,0.2,4,1.8,1.8))
  }else{
    layout(matrix(data=c(10,9,0,0,0,0,0,
                         8,7,0,0,0,0,0,
                         6,5,14,14,14,14,0,
                         4,3,0,0,0,0,0,
                         2,1,11,12,13,15,16,
                         0,17,0,0,0,0,0,
                         0,18,0,0,0,0,0), nrow=7,ncol=7,byrow=TRUE),widths=c(5*(ncol(sub.normal)/ncol(sub.tumor)),5,0.2,0.4,2,2.8,2.8),heights=c(0.2,0.2,0.2,0.2,4,1.8,1.8))
  }
  
  
  heatmap.main(cluster_T,col=jet.colors(255),margin=c(1,1,0.5,1),cexRow = 1.5, cexCol = 2,
               zlim=c(0,1),nonlab=T)
  heatmap.main(cluster_N,col=jet.colors(255),margin=c(1,1,0.5,0.5),cexRow = 1.5, cexCol = 2,nonlab=T,zlim=c(0,1))
  MultiSide.bars(data=list(lower.bound=lower.bound, cluster=factor(cutree(cluster_T$hcc,k=5)), Batch=factor(sub.Info$TCGA.BATCH),Purity=sub.Info$ABSOLUTE.purity),
                 side="colside",order=cluster_T$colInd, margins=c(0.2,1,0.2,1),
                 col=list(cluster=c("grey","white","black","pink","purple"),
                          Batch=1:length(unique(sub.Info$TCGA.BATCH)),
                          Purity=GreyCol(255),
                          lower.bound = jet.colors(255)),
                 zlim=list(lower.bound=c(0,1),Purity=c(0,1)))
  side.bars(as.matrix(sub.buffy),side="rowside",order=cluster_T$rowInd,margins=c(1,0.5,0.5,0.2),col=jet.colors(255),zlim=c(0,1),cex=1)
  side.bars2(as.matrix(sub.CGI,ncol=1),side="rowside",order=cluster_T$rowInd,margins=c(1,0.5,0.5,0.2),cex=1)
  side.bars(sub.segmentation,side="rowside",order=cluster_T$rowInd,margins=c(1,0.5,0.5,0.2),col=c("lightcyan","chocolate","yellow","purple","darkgreen","green","darkorange"),cex=1.5)
  AddLengend (Leve,cols=c("lightcyan","chocolate","yellow","purple","darkgreen","green","darkorange"),margins=c(3,1,0.5,1),cexlab=2,lab.las=1)
  side.bars(as.matrix(sub.REMC.WGBS),side="rowside",order=cluster_T$rowInd,margins=c(1,0.5,0.5,0.2),col=jet.colors(255),zlim=c(0,1),cex=1.5)
  if(paste(CT,collapse="_") %in% "2") {
    side.bars(matrix(sub.pam,ncol=1),side="colside",order=cluster_T$colInd,margins=c(0.2,1,0.2,1),col=c(1:length(unique(sub.Info$pam50))),cex=1)
    AddLengend (c("missing",levels(sub.Info$pam50)),cols=c(1:length(unique(sub.pam))),margins=c(3,1,0.5,1),cexlab=2,lab.las=1)
  }else if( length(CT) > 1){
    side.bars(matrix(sub.pam,ncol=1),side="colside",order=cluster_T$colInd,margins=c(0.2,1,0.2,1),col=c(1:length(unique(sub.pam))),cex=1)
    AddLengend (c(paste0(CT,unique(sub.pam))),cols=c(1:length(unique(sub.pam))),margins=c(3,1,0.5,1),cexlab=2,lab.las=1)
  }
  side.bars2(as.matrix(sub.TF.matrix),side="rowside",order=cluster_T$rowInd,margins=c(1,0.5,0.5,0.2),cex=1.5)
  side.bars(as.matrix(Mean.TF.sample),side="colside",order=cluster_T$colInd,margins=c(0.2,1,0.2,1),col=jet.colors(255),zlim=c(0,1),cex=1.5)
  side.bars(as.matrix(sub.ImpInfo),side="colside",order=cluster_T$colInd,margins=c(0.2,1,0.2,1),col=c("lightcyan","white","black"),cex=1.5)
  dev.off()
  rm(cluster_T,cluster_N,sub.segmentation,sub.TF.matrix,Mean.TF.sample,sub.ImpInfo,sub.CGI,sub.REMC.WGBS,TF.cluster,sub.buffy,tmp,sub.probes)
  gc()
}




