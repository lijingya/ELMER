context("Get TF")

test_that("Correclty shows TF if top5 TFs cotinas any member of the motif TF family", {
  data(elmer.data.example, envir = environment())
  enriched.motif <- list("P53_HUMAN.H10MO.B" = c("cg00329272", "cg10097755", 
                                                 "cg08928189", "cg17153775",
                                                 "cg21156590", "cg19749688", 
                                                 "cg12590404", "cg24517858", 
                                                 "cg00329272", "cg09010107", 
                                                 "cg15386853", "cg10097755", 
                                                 "cg09247779", "cg09181054"))
  suppressMessages({
    TF <- get.TFs(data, 
                  enriched.motif, 
                  TFs=data.frame(external_gene_name=c("TP53",
                                                      "TP63",
                                                      "TP73",
                                                      "DLX6",
                                                      "DMRT1"
                  ),
                  ensembl_gene_id= c("ENSG00000141510",
                                     "ENSG00000073282",
                                     "ENSG00000078900",
                                     "ENSG00000006377",
                                     "ENSG00000137090"),
                  stringsAsFactors = FALSE),
                  group.col = "shortLetterCode",
                  group1 = "TP",
                  group2 = "NT",
                  label="hypo")
  })

  tf.family <- createMotifRelevantTfs()  
  expect_true(TF$potential.TF.family %in% tf.family$P53_HUMAN.H10MO.B)
  expect_true(TF$top.potential.TF.family %in% TF$top_5percent)
  expect_true(TF$top.potential.TF.family %in% TF$potential.TF.family)
  expect_true(is.na(TF$top.potential.TF.subfamily))
  expect_true(TF$potential.TF.family %in% TF$top_5percent)
  
})  

test_that("Shows NA if top5 TFs does not include any member of the motif TF family", {
  data(elmer.data.example, envir = environment())
  enriched.motif <- list("P53_HUMAN.H10MO.B" = c("cg00329272", "cg10097755", 
                                                 "cg08928189", "cg17153775",
                                                 "cg21156590", "cg19749688", 
                                                 "cg12590404", "cg24517858", 
                                                 "cg00329272", "cg09010107", 
                                                 "cg15386853", "cg10097755", 
                                                 "cg09247779", "cg09181054"))
  suppressMessages({
    TF <- get.TFs(data, enriched.motif, label = "hypo",
                  group.col = "shortLetterCode",
                  group1 = "TP",
                  group2 = "NT")  
  })

  tf.family <- createMotifRelevantTfs()  
  human.tf <- getTF()
  # Check if top5 has 5% elements that TF from the object
  expect_equal(floor(sum(human.tf$ensembl_gene_id %in% rownames(getExp(data))) * 0.05), 
               length(unlist(strsplit(as.character(TF$top_5percent),";"))))
  if(!TF$top_5percent %in% tf.family$P53_HUMAN.H10MO.B){
    expect_true(is.na(TF$top.potential.TF.family))
    expect_true(is.na(TF$potential.TF.family))
    expect_true(is.na(TF$top.potential.TF.subfamily))
    expect_true(is.na(TF$potential.TF.subfamily))
  }
})  

test_that("Test if the results is right", {
  
  # We will create the data where whe have 3 cases: 
  # 1) no changes in expression
  # 2) Unmethylated group has lower TF expression
  # 3) Unmethylated group has a higher TF expression
  # 4) Unmethylated group has highest TF expression
  # The case 4 is the potential TF
  
  # We have the 3 cases for 6 patients
  exp <- t(data.frame("ENSG00000141510" = c(1,1,1,1,1,1), # No change in expression
                      "ENSG00000073282" = c(0,0,0,1,1,1), # Change in the other direction
                      "ENSG00000135776" = c(0.2,0.4,0.6,0.8,0.9,1), # raw p-value should be higher than the best case
                      "ENSG00000078900" = c(1,1,1,0,0,0))) # Should be true
  colnames(exp) <- c(as.character(1:6))
  exp <- makeSummarizedExperimentFromGeneMatrix(exp, genome = "hg19")
  
  # First 3 patients are Unmethylated
  met <- t(data.frame("cg00329272" = c(0,0,0,1,1,1)))
  colnames(met) <- c(as.character(1:6))
  met <- makeSummarizedExperimentFromDNAMethylation(met, met.platform = "450k", genome = "hg19")  
  
  colData <- data.frame(sample = as.character(1:6), 
                        group = c(rep("g2",3),rep("g1",3)),
                        row.names =  as.character(1:6))
  # Create datas
  data <- createMAE(exp,met, genome = "hg19", colData =  colData)
  
  enriched.motif <- list("P53_HUMAN.H10MO.B" = c("cg00329272"))
  

  suppressMessages({
    TF <- get.TFs(data, 
                  enriched.motif, 
                  group.col = "group",
                  group1 = "g1",
                  group2 = "g2",  
                  TFs = data.frame(external_gene_name=c("TP53", "TP63","TP73","ABCB10"),
                                   ensembl_gene_id= c("ENSG00000141510",
                                                      "ENSG00000073282",
                                                      "ENSG00000078900",
                                                      "ENSG00000135776"),
                                   stringsAsFactors = FALSE),
                  label = "hypo")
  })

  
  expect_true(TF$potential.TF.family == "TP73")
  expect_true(TF$top.potential.TF.family == "TP73")
  expect_true(TF$top_5percent == "TP73")
  
  # Changing percentage to 50% (split in half: 3 samples as methylated and 3 as unmethylated)
  # Will give us the same result

  suppressMessages({
    TF <- get.TFs(data, 
                  enriched.motif, 
                  minSubgroupFrac = 0.5, 
                  group.col = "group",
                  group1 = "g1",
                  group2 = "g2",
                  TFs = data.frame(external_gene_name=c("TP53", "TP63","TP73","ABCB10"),
                                   ensembl_gene_id= c("ENSG00000141510",
                                                      "ENSG00000073282",
                                                      "ENSG00000078900",
                                                      "ENSG00000135776"),
                                   stringsAsFactors = FALSE),
                  label = "hypo")
  })

  expect_true(TF$potential.TF.family == "TP73")
  expect_true(TF$top.potential.TF.family == "TP73")
  expect_true(TF$top_5percent == "TP73")
  
  
  # Changing the order should give the right gene
  exp <- t(data.frame("ENSG00000078900" = c(1,1,1,1,1,1), # No change in expression
                      "ENSG00000073282" = c(0,0,0,1,1,1), # Change in the other direction
                      "ENSG00000135776" = c(0.2,0.4,0.6,0.8,0.9,1), # raw p-value should be higher than the best case
                      "ENSG00000141510" = c(1,1,1,0,0,0))) # Should be true
  colnames(exp) <- c(as.character(1:6))
  exp <- makeSummarizedExperimentFromGeneMatrix(exp, genome = "hg19")
  
  # First 3 patients are Unmethylated
  met <- t(data.frame("cg00329272" = c(0,0,0,1,1,1)))
  colnames(met) <- c(as.character(1:6))
  met <- makeSummarizedExperimentFromDNAMethylation(met, met.platform = "450k", genome = "hg19")  
  
  colData <- data.frame(sample = as.character(1:6), 
                        group = c(rep("g2",3),rep("g1",3)),
                        row.names =  as.character(1:6))
  # Create datas
  data <- createMAE(exp,met, genome = "hg19", colData = colData)
  
  enriched.motif <- list("P53_HUMAN.H10MO.B" = c("cg00329272"))

  suppressMessages({
    TF <- get.TFs(data, 
                  enriched.motif, 
                  group.col = "group",
                  group1 = "g1",
                  group2 = "g2",
                  TFs = data.frame(external_gene_name=c("TP53", "TP63","TP73","ABCB10"),
                                   ensembl_gene_id= c("ENSG00000141510",
                                                      "ENSG00000073282",
                                                      "ENSG00000078900",
                                                      "ENSG00000135776"),
                                   stringsAsFactors = FALSE),
                  label = "hypo")
  })

  expect_true(TF$potential.TF.family == "TP53")
  expect_true(TF$top.potential.TF.family == "TP53")
  expect_true(TF$top_5percent == "TP53")
  
  # Changing percentage to 50% (split in half: 3 samples as methylated and 3 as unmethylated)
  # Will give us the same result

  suppressMessages({
    TF <- get.TFs(data, 
                  enriched.motif, 
                  minSubgroupFrac = 0.5, 
                  group.col = "group",
                  group1 = "g1",
                  group2 = "g2",
                  TFs = data.frame(external_gene_name=c("TP53", "TP63","TP73","ABCB10"),
                                   ensembl_gene_id= c("ENSG00000141510",
                                                      "ENSG00000073282",
                                                      "ENSG00000078900",
                                                      "ENSG00000135776"),
                                   stringsAsFactors = FALSE),
                  label = "hypo")
  })
  expect_true(TF$potential.TF.family == "TP53")
  expect_true(TF$top.potential.TF.family == "TP53")
  expect_true(TF$top_5percent == "TP53")
})

test_that("It creates a PDF with the TF ranking plot", {
  expect_true(file.exists("TFrankPlot_family/P53_HUMAN.H10MO.B.TFrankPlot.pdf"))
  expect_true(file.exists("TFrankPlot_subfamily/P53_HUMAN.H10MO.B.TFrankPlot.pdf"))
  unlink("TFrankPlot_family",recursive = TRUE, force = TRUE)
  unlink("TFrankPlot_subfamily",recursive = TRUE, force = TRUE)
  unlink("getTF.hypo.significant.TFs.with.motif.summary.csv",recursive = TRUE, force = TRUE)
  unlink("subfamily.motif.relevant.TFs.rda",recursive = TRUE, force = TRUE)
  unlink("HumanTF.rda",recursive = TRUE, force = TRUE)
  unlink("family.motif.relevant.TFs.rda",recursive = TRUE, force = TRUE)
  unlink("getTF.hypo.TFs.with.motif.pvalue.rda",recursive = TRUE, force = TRUE)
})
  