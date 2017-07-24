context("Testing get.diff.meth")

test_that("The directions should change if we change the groups", {
  data <- ELMER:::getdata("elmer.data.example")
  Hypo.probe.1 <- get.diff.meth(data, 
                                minSubgroupFrac = 1,
                                diff.dir="hypo",
                                group.col = "definition", 
                                group1 = "Primary solid Tumor", 
                                group2 = "Solid Tissue Normal",
                                sig.dif = 0.1) # get hypomethylated probes
  Hyper.probe.1 <- get.diff.meth(data, 
                                 minSubgroupFrac = 1,
                                 diff.dir="hyper",
                                 group.col = "definition", 
                                 group1 = "Solid Tissue Normal", 
                                 group2 = "Primary solid Tumor",
                                 sig.dif = 0.1) # get hypomethylated probes
  expect_equal(Hyper.probe.1$probe,Hypo.probe.1$probe)
  expect_true(min(Hypo.probe.1$adjust.p) >= 0)
  expect_true(max(Hypo.probe.1$adjust.p) <= 1)
  expect_true(min(Hypo.probe.1$pvalue) >= 0)
  expect_true(max(Hypo.probe.1$pvalue) <= 1)
  
  expect_equal(round(Hyper.probe.1[,3],digits = 2),-round(Hypo.probe.1[,3],digits = 2))
  
  mean.tp <- rowMeans(assay(getMet(data)[Hyper.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Primary solid Tumor"]))
  mean.nt <- rowMeans(assay(getMet(data)[Hyper.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Solid Tissue Normal"]))
  expect_equal(Hyper.probe.1[1,3] > 0, (mean.nt - mean.tp)[[1]] > 0)
  
  mean.tp <- rowMeans(assay(getMet(data)[Hypo.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Primary solid Tumor"]))
  mean.nt <- rowMeans(assay(getMet(data)[Hypo.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Solid Tissue Normal"]))
  expect_equal(Hypo.probe.1[1,3] < 0, (mean.tp - mean.nt )[[1]] < 0)
  
  Hyper.probe.2 <- get.diff.meth(data, 
                                 minSubgroupFrac = 1,
                                 diff.dir="hyper",
                                 group.col = "definition", 
                                 group1 = "Primary solid Tumor", 
                                 group2 = "Solid Tissue Normal",
                                 sig.dif = 0.1) # get hypomethylated probes
  Hypo.probe.2 <- get.diff.meth(data, 
                                minSubgroupFrac = 1,
                                diff.dir="hypo",
                                group.col = "definition", 
                                group1 = "Solid Tissue Normal", 
                                group2 = "Primary solid Tumor",
                                sig.dif = 0.1) # get hypomethylated probes
  expect_equal(Hyper.probe.2$probe,Hypo.probe.2$probe)
  expect_equal(round(Hyper.probe.2[,3],digits = 2),-round(Hypo.probe.2[,3],digits = 2))
  
  mean.tp <- rowMeans(assay(getMet(data)[Hyper.probe.2$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Primary solid Tumor"]))
  mean.nt <- rowMeans(assay(getMet(data)[Hyper.probe.2$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Solid Tissue Normal"]))
  expect_equal(Hyper.probe.2[1,3] > 0, (mean.tp - mean.nt)[[1]] > 0)
  
  mean.tp <- rowMeans(assay(getMet(data)[Hypo.probe.2$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Primary solid Tumor"]))
  mean.nt <- rowMeans(assay(getMet(data)[Hypo.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Solid Tissue Normal"]))
  expect_equal(Hypo.probe.2[1,3] < 0, (mean.nt-mean.tp )[[1]] < 0)
  
})

test_that("The test argument can be changed", {
  data <- ELMER:::getdata("elmer.data.example")

  suppressMessages({
    Hypo.probe.1 <- get.diff.meth(data, 
                                  minSubgroupFrac = 1,
                                  diff.dir="hypo",
                                  test = t.test,   
                                  group.col = "definition",
                                  group1 = "Primary solid Tumor", 
                                  group2 = "Solid Tissue Normal",
                                  sig.dif = 0.1) # get hypomethylated probes
    Hyper.probe.1 <- get.diff.meth(data, 
                                   minSubgroupFrac = 1,
                                   diff.dir="hyper",
                                   test = t.test,   
                                   group.col = "definition", 
                                   group1 = "Solid Tissue Normal", 
                                   group2 = "Primary solid Tumor",
                                   sig.dif = 0.1) # get hypomethylated probes
  })
  expect_equal(Hyper.probe.1$probe,Hypo.probe.1$probe)
  expect_equal(round(Hyper.probe.1[,3],digits = 2),-round(Hypo.probe.1[,3],digits = 2))
  
  mean.tp <- rowMeans(assay(getMet(data)[Hyper.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Primary solid Tumor"]))
  mean.nt <- rowMeans(assay(getMet(data)[Hyper.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Solid Tissue Normal"]))
  expect_equal(Hyper.probe.1[1,3] > 0, (mean.nt - mean.tp)[[1]] > 0)
  
  mean.tp <- rowMeans(assay(getMet(data)[Hypo.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Primary solid Tumor"]))
  mean.nt <- rowMeans(assay(getMet(data)[Hypo.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Solid Tissue Normal"]))
  expect_equal(Hypo.probe.1[1,3] < 0, (mean.tp - mean.nt )[[1]] < 0)
  

  suppressMessages({
    Hyper.probe.2 <- get.diff.meth(data, 
                                   minSubgroupFrac = 1,
                                   diff.dir="hyper",
                                   group.col = "definition", 
                                   test = t.test,   
                                   group1 = "Primary solid Tumor", 
                                   group2 = "Solid Tissue Normal",
                                   sig.dif = 0.1) # get hypomethylated probes
    Hypo.probe.2 <- get.diff.meth(data, 
                                  minSubgroupFrac = 1,
                                  diff.dir="hypo",
                                  test = t.test,   
                                  group.col = "definition", 
                                  group1 = "Solid Tissue Normal", 
                                  group2 = "Primary solid Tumor",
                                  sig.dif = 0.1) # get hypomethylated probes
  })
  expect_equal(Hyper.probe.2$probe,Hypo.probe.2$probe)
  expect_equal(round(Hyper.probe.2[,3],digits = 2),-round(Hypo.probe.2[,3],digits = 2))
  
  mean.tp <- rowMeans(assay(getMet(data)[Hyper.probe.2$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Primary solid Tumor"]))
  mean.nt <- rowMeans(assay(getMet(data)[Hyper.probe.2$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Solid Tissue Normal"]))
  expect_equal(Hyper.probe.2[1,3] > 0, (mean.tp - mean.nt)[[1]] > 0)
  
  mean.tp <- rowMeans(assay(getMet(data)[Hypo.probe.2$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Primary solid Tumor"]))
  mean.nt <- rowMeans(assay(getMet(data)[Hypo.probe.1$probe[1],colData(data)[sampleMap(data)[sampleMap(data)$assay == "DNA methylation","primary"],"definition"] == "Solid Tissue Normal"]))
  expect_equal(Hypo.probe.2[1,3] < 0, (mean.nt-mean.tp )[[1]] < 0)
})

test_that("It threats correclty NAs and thrseholds", {
  data <- ELMER:::getdata("elmer.data.example")

  suppressMessages({
    diff.probes <- get.diff.meth(data, 
                                 minSubgroupFrac = 1,
                                 diff.dir="hypo",
                                 test = t.test,   
                                 group.col = "definition",
                                 group1 = "Primary solid Tumor", 
                                 group2 = "Solid Tissue Normal",
                                 sig.dif = 0.1,
                                 pvalue = 0.0) # get hypomethylated probes
    expect_true(nrow(diff.probes) == 0)
    diff.probes <- get.diff.meth(data, 
                                 minSubgroupFrac = 1,
                                 diff.dir="hypo",
                                 test = t.test,   
                                 group.col = "definition",
                                 group1 = "Primary solid Tumor", 
                                 group2 = "Solid Tissue Normal",
                                 sig.dif = 0.0,
                                 pvalue = 1.1) # get hypomethylated probes
  })
  expect_equal(nrow(diff.probes), nrow(diff.probes))
})

