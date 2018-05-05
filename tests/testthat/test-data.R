context("Checking if data imported is correctly created")

test_that("Family of TF is correctly created from HOCOMOCO (TFClass database)", {
  # Create list of TF in the same family
  family <- createMotifRelevantTfs()
  # TP53 TP63 ad TP73 are in the same family
  expect_true(all(c("TP53", "TP63", "TP73") %in% family$P53_HUMAN.H11MO.0.A))
  expect_true(all(c("TP53", "TP63", "TP73") %in% family$P63_HUMAN.H11MO.0.A))
  expect_true(all(c("TP53", "TP63", "TP73") %in% family$P73_HUMAN.H11MO.0.A))
  expect_equal(length(family$P53_HUMAN.H11MO.0.A), 3)
  expect_equal(length(family$P63_HUMAN.H11MO.0.A), 3)
  expect_equal(length(family$P73_HUMAN.H11MO.0.A), 3)
})

test_that("Get list of human TF from ELMER.data", {
  tf <- getTF()
  expect_true(all(c("TP53", "TP63", "TP73") %in% tf$external_gene_name))
})

test_that("Mapping from entrez gene ID to emsemble gene ID is right", {
  genes <- c("100887754","100873766","100874231")
  mapping <- get.GRCh(genome = "hg19",genes)
  expect_equal(mapping[mapping$entrezgene == genes[1],]$ensembl_gene_id,"ENSG00000243300" )
  expect_equal(mapping[mapping$entrezgene == genes[2],]$ensembl_gene_id,"ENSG00000252952" )
  expect_equal(mapping[mapping$entrezgene == genes[3],]$ensembl_gene_id,"ENSG00000227213" )
})