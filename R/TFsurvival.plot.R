#' Creates survival plot of based on the expression of a TF
#' @description 
#' This function will create a survival plot for the samples with higher, midium, low expression
#' of a given transcription factor. 
#' By defau;t samples with higher expression are the top 30% and the lower expression the bottom 30%.
#' @param data A multi assay Experiment with clinical data in the  phenotypic data matrix
#'  containing the following columns: vital_status, days_to_last_follow_up and days_to_death. Default from GDC and TCGAbiolinks
#' @param TF A gene symbol
#' @param xlim Limit x axis showed in plot
#' @param percentage   A number ranges from 0 to 1 specifying the percentage of samples in the 
#' higher and lower expression groups. Default is 0.3
#' @param save Save plot as PDF
#' @importFrom TCGAbiolinks TCGAanalyze_survival
#' @export
TFsurvival.plot <- function(data,
                            TF,
                            xlim = NULL,
                            percentage = 0.3,
                            save = TRUE){
  if(!all(c("vital_status", "days_to_last_follow_up","days_to_death") %in% colnames(colData(data)))){
    message("colData must have the following columns: vital_status,days_to_last_follow_up, days_to_death")
    return(NULL)
  }
  # For the transcription factor, gets it getGeneID
  gene <- getGeneID(data,symbol=TF)
  # Get the expression values for the genes.
  # (getExp is a ELMER function)
  exp <- as.vector(assay(getExp(data)[gene,]))
  names(exp) <- colnames(getExp(data))
  exp <- sort(exp)
  
  # Get the names of the 30% patients with lower expression
  lower <- names(head(exp, n = ceiling(length(exp) * percentage)))
  
  # Get the names of the 30% patients with higher expression
  higher <- names(tail(exp, n = ceiling(length(exp) * percentage)))
  
  df <- colData(data)
  # Create the labels for each sample
  df$tf_groups <- "medium"
  low.idx <-  sampleMap(data)[sampleMap(data)$colname %in% lower,"primary"]
  df[low.idx,"tf_groups"] <- "low"
  high.idx <-  sampleMap(data)[sampleMap(data)$colname %in% higher,"primary"]
  df[high.idx,"tf_groups"] <- "high"
  
  filename <-  NULL
  if(save) filename <- paste0(TF,"_survival.pdf")
  # Use TCGAbiolinks to create the survival curve
  TCGAanalyze_survival(df,
                       "tf_groups",
                       legend=paste0(TF," Exp level"),
                       filename = filename,
                       xlim = xlim,
                       conf.int = FALSE,
                       risk.table = FALSE)
}
