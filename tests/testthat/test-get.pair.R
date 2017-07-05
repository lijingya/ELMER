context("Checking get pair function")

test_that("Function uses correctly the permu.dir", {
  data(elmer.data.example, envir = environment())
  sink("/dev/null");
  
  nearGenes <- GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
                            geneAnnot=getExp(data))
  Hypo.pair <- get.pair(data = data,
                        group.col = "definition", 
                        group1 = "Primary solid Tumor", 
                        group2 = "Solid Tissue Normal",
                        nearGenes = nearGenes,
                        permu.size = 5,
                        pvalue =  0.001,
                        Pe = 0.2,
                        dir.out="./",
                        permu.dir = "permu_test",
                        label = "hypo")
  # Folder was crreated correcly
  expect_true(file.exists("permu_test/permu.rda"))
  expect_true(ncol(get(load("permu_test/permu.rda"))) == 5)
  # Result correctly use a gene from the nearGene list
  expect_true(all(Hypo.pair[Hypo.pair$Probe %in% "cg00329272" ,]$GeneID %in% nearGenes$cg00329272$GeneID))
  expect_true(all(Hypo.pair$Pe < 0.2))
  expect_true(min(Hypo.pair$Pe) >= 0)
  expect_true(max(Hypo.pair$Pe) <= 1)
  expect_true(min(Hypo.pair$Raw.p) >= 0)
  expect_true(max(Hypo.pair$Raw.p) <= 1)
    
  # If we add one more probe the value should be saved
  Hypo.pair <- get.pair(data=data,
                        nearGenes=nearGenes,
                        permu.size = 6,
                        group.col = "definition", 
                        group1 = "Primary solid Tumor", 
                        group2 = "Solid Tissue Normal",
                        pvalue =  0.2,
                        dir.out="./",
                        permu.dir = "permu_test",
                        label= "hypo")
  # Folder was crreated correcly
  expect_true(file.exists("permu_test/permu.rda"))
  expect_true(ncol(get(load("permu_test/permu.rda"))) == 6)
  
  # If we reduce the number of probes
  Hypo.pair <- get.pair(data=data,
                        nearGenes=nearGenes,
                        group.col = "definition", 
                        group1 = "Primary solid Tumor", 
                        group2 = "Solid Tissue Normal",
                        permu.size=5,
                        pvalue =  0.02,
                        dir.out="./",
                        label= "hypo")
  # Pvalue filter is working
  expect_true(nrow(Hypo.pair) == 0)
  
  # If we add new genes
  nearGenes <- GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755","cg22396959"),],
                            geneAnnot=getExp(data))
  Hypo.pair <- get.pair(data=data,
                        group.col = "definition", 
                        group1 = "Primary solid Tumor", 
                        group2 = "Solid Tissue Normal",
                        nearGenes=nearGenes,
                        permu.size=7,
                        pvalue =  0.2,
                        dir.out="./",
                        permu.dir = "permu_test",
                        label= "hypo")
})

test_that("Gene expression is calculated", {
  data(elmer.data.example, envir = environment())
  sink("/dev/null");
  
  nearGenes <- GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
                            geneAnnot=getExp(data))
  Hypo.pair <- get.pair(data = data,
                        nearGenes = nearGenes,
                        permu.size = 5,
                        pvalue =  0.001,
                        Pe = 0.2,
                        diffExp = TRUE,
                        group.col = "definition", 
                        group1 = "Primary solid Tumor", 
                        group2 = "Solid Tissue Normal",
                        dir.out="./",
                        permu.dir = "permu_test",
                        label = "hypo")
  expect_true(any(grepl("log2FC", colnames(Hypo.pair))))
  expect_true(any(grepl("pvalue", colnames(Hypo.pair))))
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
test_that("Ramdom probe sleection is the same for every run", {
  probes <- paste0("cg",000000:450000)
  set.seed(200); probes.permu <- sample(probes, size = 10000, replace = FALSE)
  set.seed(200); probes.permu.rep <- sample(probes, size = 10000, replace = FALSE)
  expect_true(all(probes.permu == probes.permu.rep))
  
  data(elmer.data.example, envir = environment())
  permu <- get.permu(data = data,
                     permu.dir = "test_permu_1",
                     geneID=rownames(getExp(data))[1],
                     rm.probes=c("cg00329272","cg10097755"),
                     permu.size=51)
  permu <- get.permu(data = data,
                     permu.dir = "test_permu_2",
                     geneID=rownames(getExp(data))[1],
                     rm.probes=c("cg00329272","cg10097755"),
                     permu.size=51)
  probes.permu <- colnames(get(load("test_permu_1/permu.rda")))
  probes.permu.rep <-  colnames(get(load("test_permu_2/permu.rda")))
  expect_true(all(probes.permu == probes.permu.rep))
  unlink("test_permu_2",recursive = TRUE, force = TRUE)
  unlink("test_permu_1",recursive = TRUE, force = TRUE)
  unlink("permu_test",recursive = TRUE, force = TRUE)
})
