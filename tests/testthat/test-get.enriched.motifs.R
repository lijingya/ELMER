context("Check get enriched motif function")

test_that("get enriched motif function returns the expected result", {
  data(elmer.data.example, envir = environment())
  bg <- rownames(getMet(data))
  probes <- bg[1:20]
  
  # In this case MAFG_HUMAN.H10MO.C is the enriched motif
  # 1) SP2_HUMAN.H10MO.C, has the motif for all probes, but as it has for all bg probes
  # it will be conisered false positve
  # 2) Has 1 for all the probes and 0 for all background (our best case)
  # 3) Has 0 for all cases (should not be selected)
  Probes.motif <- data.frame("SP2_HUMAN.H10MO.C" = rep(1, length(bg)),
                             "MAFG_HUMAN.H10MO.C" =  c(rep(1, length(bg)/2),rep(0, length(bg)/2 + 1)),
                             "NR2E1_HUMAN.H10MO.D" = rep(0, length(bg)),
                             "TBX15_HUMAN.H10MO.D" = rep(0, length(bg)))
  rownames(Probes.motif) <- bg
  Probes.motif[probes,4] <- 1
  Probes.motif[1,] <- c(0,0,1,0) # The case before will give an execption

  suppressMessages({
    enriched.motif <- get.enriched.motif(probes.motif=Probes.motif,
                                         probes=probes,lower.OR = 0.1,
                                         min.motif.quality = "D",
                                         background.probes = bg,
                                         label="hypo")
  })
  # In this case MAFG_HUMAN.H10MO.C is the enriched motif
  # 1) SP2_HUMAN.H10MO.C, has the motif for all probes, but as it has for all bg probes
  # it will be conisered false positve
  # 2) Has 1 for all the probes and 0 for all background (our best case)
  # 3) Has 0 for all cases (should not be selected)
  Probes.motif <- data.frame("SP2_HUMAN.H10MO.C" = rep(1, length(bg)),
                             "MAFG_HUMAN.H10MO.C" =  c(rep(1, length(bg)/2),rep(0, length(bg)/2 + 1)),
                             "NR2E1_HUMAN.H10MO.D" = rep(0, length(bg)))
  rownames(Probes.motif) <- bg
  Probes.motif[1,] <- c(0,0,1) # The case before will give an execption
  

  suppressMessages({
    enriched.motif <- get.enriched.motif(probes.motif=Probes.motif,
                                         probes=probes,
                                         min.motif.quality = "D",
                                         background.probes = bg,
                                         label="hypo")
  })
  expect_equal(names(enriched.motif), "MAFG_HUMAN.H10MO.C")
  expect_true(all(enriched.motif[[1]] %in% probes))
})

test_that("min.incidence works", {
  data(elmer.data.example, envir = environment())
  bg <- rownames(getMet(data))
  probes <- bg[1:20]
  
  # In this case MAFG_HUMAN.H10MO.C is the enriched motif
  # 1) SP2_HUMAN.H10MO.C, has the motif for all probes, but as it has for all bg probes
  # it will be conisered false positve
  # 2) Has 1 for all the probes and 0 for all background (our best case)
  # 3) Has 0 for all cases (should not be selected)
  Probes.motif <- data.frame("SP2_HUMAN.H10MO.C" = rep(1, length(bg)),
                             "MAFG_HUMAN.H10MO.C" =  c(rep(1, length(bg)/2),rep(0, length(bg)/2 + 1)),
                             "NR2E1_HUMAN.H10MO.D" = rep(0, length(bg)),
                             "TBX15_HUMAN.H10MO.D" = rep(0, length(bg)))
  rownames(Probes.motif) <- bg
  Probes.motif[probes,4] <- 1
  Probes.motif[1,] <- c(0,0,1,0) # The case before will give an execption

  suppressMessages({
    
    enriched.motif <- get.enriched.motif(probes.motif=Probes.motif,
                                         probes=probes,
                                         min.incidence = 20,
                                         lower.OR = 0.1,
                                         min.motif.quality = "D",
                                         background.probes = bg,
                                         label="hypo")
  })
  expect_true(length(enriched.motif) == 0)

  suppressMessages({
    
    enriched.motif <- get.enriched.motif(probes.motif=Probes.motif,
                                         probes=probes,
                                         min.incidence = 0,
                                         lower.OR = 0.0,
                                         min.motif.quality = "D",
                                         background.probes = bg,
                                         label="hypo")
  })
  expect_true(length(enriched.motif) == ncol(Probes.motif))

  suppressMessages({
    enriched.motif <- get.enriched.motif(probes.motif=Probes.motif,
                                         probes=probes,
                                         min.incidence = 0,
                                         lower.OR = 0.0,
                                         min.motif.quality = "B",
                                         background.probes = bg,
                                         label="hypo")
  })
  expect_true(length(enriched.motif) == 0)
  unlink("hypo.quality*",recursive = TRUE, force = TRUE)
  unlink("hypo.all*",recursive = TRUE, force = TRUE)
  unlink("getMotif.hypo.enriched.motifs.rda",recursive = TRUE, force = TRUE)
  unlink("getMotif.hypo.motif.enrichment.csv",recursive = TRUE, force = TRUE)
})