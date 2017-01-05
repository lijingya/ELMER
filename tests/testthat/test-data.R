context("Checking if data imported is correctly created")

test_that("Family of TF is correctly created from HOCOMOCO (TFClass database)", {
  # Create list of TF in the same family
  family <- createMotifRelevantTfs()
  # TP53 TP63 ad TP73 are in the same family
  expect_true(all(c("TP53", "TP63", "TP73") %in% family$P53_HUMAN.H10MO.B))
  expect_true(all(c("TP53", "TP63", "TP73") %in% family$P63_HUMAN.H10MO.A))
  expect_true(all(c("TP53", "TP63", "TP73") %in% family$P73_HUMAN.H10MO.A))
  expect_equal(length(family$P53_HUMAN.H10MO.B), 3)
  expect_equal(length(family$P63_HUMAN.H10MO.A), 3)
  expect_equal(length(family$P73_HUMAN.H10MO.A), 3)
})

test_that("Get list of TF from uniprot database", {
  tf <- getTF()
  expect_true(all(c("TP53", "TP63", "TP73") %in% tf$external_gene_name))
  expect_true(all(c("ensembl_gene_id", "external_gene_name") %in% colnames(tf)))
})
