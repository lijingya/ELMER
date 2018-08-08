context("Checking getRandomPair function")
test_that("Links are as expected", {
  
  data <- ELMER:::getdata("elmer.data.example")
  nearGenes <- GetNearGenes(TRange=getMet(data),
                            geneAnnot=getExp(data))
  links <- rbindlist(nearGenes)
  links <-  links[sample(1:nrow(links),250)] # get 250 random links
  random.pairs <- getRandomPairs(as.data.frame(links))
  
  library(data.table)
  DT <- data.table(links[order(links$Side),], key="Target")
  sig.pairs.links <- DT[, list(col1=paste(Side,collapse = ",")),by=Target]
  
  DT <- data.table(random.pairs[order(random.pairs$Side),], key="Probe")
  random.pairs.links <- DT[, list(col1=paste(Side,collapse = ",")),by=Probe]
  
  # Same nb of probes ? 
  expect_true(length(unique(links$Target)) == length(unique(random.pairs$Probe)))
  
  # Same number of position links
  expect_true(all(table(links$Side) == table(random.pairs$Side)))
  
  # same links per probe
  expect_true(all(table(random.pairs.links$col1) == table(sig.pairs.links$col1)))
  
})
