context("Get TF")

test_that("Correclty shows TF if top5 TFs cotinas any member of the motif TF family", {
  data(elmer.data.example)
  enriched.motif <- list("P53_HUMAN.H10MO.B" = c("cg00329272", "cg10097755", 
                                                 "cg08928189", "cg17153775",
                                                 "cg21156590", "cg19749688", 
                                                 "cg12590404", "cg24517858", 
                                                 "cg00329272", "cg09010107", 
                                                 "cg15386853", "cg10097755", 
                                                 "cg09247779", "cg09181054"))
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
                label="hypo")
  expect_true(TF$potential.TFs %in% tf.family$P53_HUMAN.H10MO.B)
  expect_true(TF$top.potential.TF %in% TF$top_5percent)
  expect_true(TF$top.potential.TF %in% TF$potential.TFs)
  expect_true(TF$potential.TFs %in% TF$top_5percent)
  
  # In this example TP53 is the one with lower raw p-value, so it shoud be the one
  # in top_5percent and in the other columns 
  expect_true(TF$potential.TFs == "TP73")
  expect_true(TF$top.potential.TF == "TP73")
  expect_true(TF$top_5percent == "TP73")
})  

test_that("Shows NA if top5 TFs does not include any member of the motif TF family", {
  data(elmer.data.example)
  enriched.motif <- list("P53_HUMAN.H10MO.B" = c("cg00329272", "cg10097755", 
                                                 "cg08928189", "cg17153775",
                                                 "cg21156590", "cg19749688", 
                                                 "cg12590404", "cg24517858", 
                                                 "cg00329272", "cg09010107", 
                                                 "cg15386853", "cg10097755", 
                                                 "cg09247779", "cg09181054"))
  TF <- get.TFs(data, 
                enriched.motif, 
                label = "hypo")  
  
  tf.family <- createMotifRelevantTfs()  
  human.tf <- getTF()
  # Check if top5 has 5% elements that TF from the object
  expect_equal(floor(sum(tf$ensembl_gene_id %in% rownames(getExp(data))) * 0.05), 
               length(unlist(strsplit(as.character(TF.all$top_5percent),";"))))
  if(!TF$top_5percent %in% tf.family$P53_HUMAN.H10MO.B){
    expect_true(is.na(TF$top.potential.TF))
    expect_true(is.na(TF$potential.TFs))
  }
})  

