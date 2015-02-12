# Custom.pipe <- function(analysis="all",wd="./",Data=NULL,...){
#   if(missing(disease)) stop("Disease should be specified.\nDisease short name (such as LAML) please check https://tcga-data.nci.nih.gov/tcga/.")
#   if(analysis == "all") analysis=c("download","distal","diffMeth","pair","motif","TF.search")
#   setwd(wd)
#   args <- list(...)
#   
#   if(is.null(Data)){
#     stop("Data need to be specified.\nEither provide downloaded data folder contain your data stored as XXX.rda.")
#   }
#   
#   ## select distal enhancer probes
#   if("distal" %in% analysis){
#     cat("########################################\nSelect distal enhancer probes\n########################################")
#     params <- args[names(args) %in% c(need to be done)]
#     get.distal.en
#   }
#   #get differential DNA methylation
#   if("diffMeth" %in% analysis){
#     cat("########################################\nGet differential DNA methylation loci\n########################################")
#     params <- args[names(args) %in% c(need to be done)]
#     Get.hypo
#     Get.hyper
#   }
#   
#   #predict pair
#   if("pair" %in% analysis){
#     cat("########################################\nPredict pairs\n########################################")
#     params <- args[names(args) %in% c(need to be done)]
#     Get.pair
#   }
#   
#   # search enriched motif
#   if("motif" %in% analysis){
#     cat("########################################\nMotif search\n########################################")
#     params <- args[names(args) %in% c(need to be done)]
#     motif.enrich
#   }
#   
#   #search responsible TFs
#   if("TF.search" %in% analysis){
#     cat("########################################\nSearch responsible TFs\n########################################")
#     params <- args[names(args) %in% c(need to be done)]
#     TF.search
#   }
# }