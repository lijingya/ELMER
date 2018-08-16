context("Checking getRandomPair function")
library(plyr)
library(dplyr)
library(data.table)

test_that("Links are as expected", {
  
  data <- ELMER:::getdata("elmer.data.example")
  links <- GetNearGenes(TRange=rowRanges(getMet(data)),
                        geneAnnot=rowRanges(getExp(data)))
  links <-  links[sample(1:nrow(links),250),] # get 250 random links
  random.pairs <- getRandomPairs(links)
  
  random.pairs %>% 
    group_by(Probe) %>%
    summarize(col1=paste(sort(Side),collapse = ",")) %>%
    data.frame() -> sig.pairs.links
  links %>% 
    group_by(Target) %>%
    summarize(col1=paste(sort(Side),collapse = ",")) %>%
    data.frame() -> random.pairs.links
  
  
  # Same nb of probes ? 
  expect_true(length(unique(links$Target)) == length(unique(random.pairs$Probe)))
  
  # Same number of position links
  expect_true(all(table(links$Side) == table(random.pairs$Side)))
  
  # same links per probe
  expect_true(all(table(random.pairs.links$col1) == table(sig.pairs.links$col1)))
})
