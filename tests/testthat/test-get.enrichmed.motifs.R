context("Check get enriched motif function")

test_that("Family of TF is correctly created from HOCOMOCO (TFClass database)", {
  data(elmer.data.example)
  bg <- rownames(getMet(data))
  probes <- bg[1:20]
  Probes.motif <- data.frame("SP2_HUMAN.H10MO.C" = rep(1, length(bg)),
                             "MAFG_HUMAN.H10MO.C" =  c(rep(1, length(bg)/2),rep(0, length(bg)/2 + 1)),
                             "NR2E1_HUMAN.H10MO.D" = rep(0, length(bg)),
                             "TBX15_HUMAN.H10MO.D" = rep(0, length(bg)))
  rownames(Probes.motif) <- bg
  Probes.motif[probes,4] <- 1
  Probes.motif[1,] <- c(0,0,1,0) # The case before will give an execption

  enriched.motif <- get.enriched.motif(probes.motif=Probes.motif,
                                       probes=probes,lower.OR = 0.01,
                                       background.probes = bg,
                                       label="hypo")
  Probes.motif <- data.frame("SP2_HUMAN.H10MO.C" = rep(1, length(bg)),
                             "MAFG_HUMAN.H10MO.C" =  c(rep(1, length(bg)/2),rep(0, length(bg)/2 + 1)),
                             "NR2E1_HUMAN.H10MO.D" = rep(0, length(bg)))
  rownames(Probes.motif) <- bg
  Probes.motif[1,] <- c(0,0,1) # The case before will give an execption
  
  enriched.motif <- get.enriched.motif(probes.motif=Probes.motif,
                                       probes=probes,
                                       background.probes = bg,
                                       label="hypo")
  
})