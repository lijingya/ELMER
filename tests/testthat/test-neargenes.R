context("Get GetNearGenes")

test_that("It maps correctly to hg38", {
  geneAnnot <- getTSS(genome = "hg38")
  probe <- GRanges(seqnames = c("chr1"), 
                   row.names =  c("cg18108049"),
                   names =  c("cg18108049"),
                   range = IRanges(start = c(16058489), end= c(16058489)))
  # chr1:16010827:16062808
  # chr1:16058489:16058489
  NearbyGenes <- GetNearGenes(numFlankingGenes = 2,
                              geneAnnot = geneAnnot,
                              TRange = probe)
  res <- NearbyGenes$cg18108049
  expect_equal(res[res$Side == "L1","GeneID"], "ENSG00000184908")
  expect_equal(res[res$Side == "R1","GeneID"], "ENSG00000185519")
  expect_equal(res[res$Side == "L1","Symbol"], "CLCNKB")
  expect_equal(res[res$Side == "R1","Symbol"], "FAM131C")
})

test_that("It maps correctly to hg19", {
  geneAnnot <- getTSS(genome = "hg19")
  probe <- GRanges(seqnames = c("chr1"), 
                   row.names =  c("cg18108049"),
                   names =  c("cg18108049"),
                   range = IRanges(start = c(16058489), end= c(16058489)))
  # chr1:16010827:16062808
  # chr1:16058489:16058489
  NearbyGenes <- GetNearGenes(numFlankingGenes = 4,
                              geneAnnot = geneAnnot,
                              TRange = probe)
  res <- NearbyGenes$cg18108049
  expect_true("RP11-288I21.1" %in% res$Symbol)
  expect_true("AL121992.1" %in% res$Symbol)
  expect_true("PLEKHM2" %in% res$Symbol)
  expect_true("SLC25A34" %in% res$Symbol)
})