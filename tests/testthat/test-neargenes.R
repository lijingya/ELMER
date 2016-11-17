context("Get GetNearGenes")

test_that("It maps correctly to hg19", {
  geneAnnot <- txs(TSS=list(upstream=0, downstream=0))
  probe <- GRanges(seqnames = c("chr1"), 
                   range = IRanges(start = c(16058489), end= c(16058489)), 
                   name = c("cg18108049"))
  # chr1:16010827:16062808
  # chr1:16058489:16058489
  NearbyGenes <- GetNearGenes(geneNum=2,geneAnnot=geneAnnot,TRange=probe)
  res <- NearbyGenes$cg18108049
  expect_equal( res[res$Side == "L1","Symbol"], "PLEKHM2")
  expect_equal( res[res$Side == "R1","Symbol"], "SLC25A34")
})