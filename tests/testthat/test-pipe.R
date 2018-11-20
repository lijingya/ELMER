context("Testing TCGA.pipe")

test_that("TCGA.pipe is working", {
  data <- ELMER:::getdata("elmer.data.example")
  TCGA.pipe(disease = "LUSC",
            data = data,
            analysis = c("diffMeth","pair", "motif","TF.search"), 
            mode = "supervised",
            group.col = "definition",
            group1 = "Primary solid Tumor", 
            group2 = "Solid Tissue Normal",
            diff.dir = c("hypo"),
            dir.out = "pipe",
            sig.dif = 0.0001,
            pvalue = 1.0,
            min.incidence = 0,
            lower.OR = 0.0)
  
})


