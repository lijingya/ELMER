context("Checking get pair function")

test_that("It maps correctly to hg19", {
  data(elmer.data.example)
  nearGenes <-GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
                           geneAnnot=getExp(data))
  Hypo.pair <-get.pair(data=data,
                       nearGenes=nearGenes,
                       permu.size=5,
                       Pe = 0.2,
                       dir.out="./",
                       label= "hypo")
  
  
})

test_that("Test calculation of Pe (empirical pvalue) from Raw-pvalue is working", {
  
  # If my raw-pvalue is smaller than for other probes my Pe should be small
  # If my raw-pvalue is higher than for other probes my Pe should be higher
  # Case 1 (ENSG00000157916):smaller
  # Case 2 (ENSG00000149527) intermediarie
  # Case 3 (ENSG00000116213) higher
  U.matrix <- data.frame("GeneID" = c("ENSG00000157916","ENSG00000149527","ENSG00000116213"),
                         "Raw.p" = c(0.001, 0.01, 0.1))
  permu <- data.frame("cg13480549" = c(0.1,0.1,0.01),
                      "cg15128801"= c(0.1,0.001,0.01),
                      "cg22396959"= c(0.1,0.001,0.01),
                      "cg13918150"= c(0.1,0.001,0.01),
                      "cg26403223"= c(0.1,0.1,0.01),
                      row.names = c("ENSG00000157916","ENSG00000149527","ENSG00000116213")
  )
  Pe <- Get.Pvalue.p(U.matrix = U.matrix, permu = permu)
  expect_true(Pe[Pe$GeneID== "ENSG00000157916","Pe"] == min(Pe$Pe))
  expect_true(Pe[Pe$GeneID== "ENSG00000149527","Pe"] < max(Pe$Pe) & Pe[Pe$GeneID== "ENSG00000149527","Pe"] > min(Pe$Pe))
  expect_true(Pe[Pe$GeneID== "ENSG00000116213","Pe"] == max(Pe$Pe))
})