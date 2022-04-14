#' @title Build report for TCGA.pipe function
#' @description Build HTML report 
#' @param title HTML report title
#' @param mae.file  Absolute path to the mae used in the analysis (.rda or .rds)
#' @param group.col Group col
#' @param group1 Group 1 
#' @param group2 Group 2 
#' @param direction direction used in the analysis
#' @param dir.out Absolute path to folder with results. dir.out used in the analysis
#' @param genome Genome of reference used in the analysis
#' @param mode mode used in the analysis
#' @param minSubgroupFrac  minSubgroupFrac used in the analysis
#' @param minMetdiff minMetdiff used in the analysis
#' @param metfdr metfdr used in the analysis
#' @param permu permu used in the analysis
#' @param rawpval rawpval used in the analysis
#' @param pe pe used in the analysis
#' @param nprobes nprobes used in the analysis
#' @param lower.OR lower.OR used in the analysis
#' @param out_file Output file name (i.e report.html)
#' @param funcivar Include funcivar analysis? 
#' @export
#' @importFrom rmarkdown render
#' @examples 
#' \dontrun{
#' render_report(
#'  group.col = "TN",
#'  group1 = "Tumor",
#'  group2 = "Normal",
#'  dir.out = "~/paper_elmer/Result/BRCA/TN_Tumor_vs_Normal/hypo/",
#'  direction = "hypo",
#'  mae.file = "~/paper_elmer/Result/BRCA/BRCA_mae_hg38.rda"
#')
#' }
render_report <- function(
  title = "Report",
  mae.file,
  group.col,
  group1,
  group2,
  direction,
  dir.out,
  genome = "hg38",
  mode = "supervised",
  minSubgroupFrac = 0.2,
  minMetdiff = 0.3,
  metfdr = 0.01,
  permu = 10000,
  rawpval = 0.01,
  pe = 0.01,
  nprobes = 10,
  lower.OR = 1.1,
  out_file = file.path(getwd(),"report.html"),
  funcivar = FALSE
) {
  if(missing(dir.out)) stop("Please, set dir.out value")
  if(missing(mae.file)) stop("Please, set mae.file value")
  if(missing(group.col)) stop("Please, set mae value")
  if(missing(group1)) stop("Please, set mae value")
  if(missing(group2)) stop("Please, set mae value")
  template <- system.file("rmd", "template.Rmd", package="ELMER")
  
  message("Compiling report")
  parameters <- list(
    title = title,
    genome = genome,
    mode = mode,
    minSubgroupFrac = minSubgroupFrac,
    minMetdiff = minMetdiff,
    metfdr = metfdr,
    permu = permu,
    rawpval = rawpval,
    pe = pe,
    nprobes = nprobes,
    lower.OR = lower.OR,
    groupCol =  group.col,
    mae.file = mae.file,
    group1 = group1,
    group2 = group2,
    direction = direction,
    dir.out = dir.out,
    funcivar = funcivar
  )
  message("Saving report: ", out_file)
  rmarkdown::render(
    template,
    intermediates_dir = dirname(out_file),
    output_file = out_file,
    params = parameters
  )
  invisible(TRUE)
}