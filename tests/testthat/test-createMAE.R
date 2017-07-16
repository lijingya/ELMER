context("Testing Multi Assay Experiment creation")

test_that("The creation of a using matrices and no TCGA data with equal colnames in DNA methylation and Gene Expression", {
  # NON TCGA example: matrices has diffetrent column names
  gene.exp <- DataFrame(sample1 = c("ENSG00000141510"=2.3,"ENSG00000171862"=5.4),
                        sample2 = c("ENSG00000141510"=1.6,"ENSG00000171862"=2.3)
  )
  dna.met <- DataFrame(sample1 = c("cg14324200"=0.5,"cg23867494"=0.1),
                       sample2 =  c("cg14324200"=0.3,"cg23867494"=0.9)
  )
  sample.info <- DataFrame(sample.type = c("Normal", "Tumor"))
  rownames(sample.info) <- colnames(gene.exp)
  sink("/dev/null");
  suppressMessages({
    mae <- createMAE(exp = gene.exp, met = dna.met, colData  = sample.info, genome = "hg38") 
  })
  expect_equal(metadata(mae)$genome,"hg38")
  expect_false(metadata(mae)$TCGA)
  expect_equal(dim(getExp(mae)),c(2,2))
  expect_equal(dim(getMet(mae)),c(2,2))
  expect_equal(assay(getMet(mae)),as.matrix(dna.met))
  expect_equal(assay(getExp(mae)),as.matrix(gene.exp))
  expect_equal(colData(mae),sample.info)
  expect_true(all(sampleMap(mae)$assay %in% c("DNA methylation","Gene expression")))
})  

test_that("The creation of a using matrices and no TCGA data with different colnames in DNA methylation and Gene Expression", {
  # NON TCGA example: matrices has diffetrent column names
  gene.exp <- DataFrame(sample1.exp = c("ENSG00000141510"=2.3,"ENSG00000171862"=5.4),
                        sample2.exp = c("ENSG00000141510"=1.6,"ENSG00000171862"=2.3)
  )
  dna.met <- DataFrame(sample1.met = c("cg14324200"=0.5,"cg23867494"=0.1),
                       sample2.met =  c("cg14324200"=0.3,"cg23867494"=0.9)
  )
  sample.info <- DataFrame(sample.type = c("Normal", "Tumor"))
  rownames(sample.info) <- c("sample1","sample2")
  sampleMap <- DataFrame(primary = c("sample1","sample1","sample2","sample2"), 
                         colname = c("sample1.exp","sample1.met","sample2.exp","sample2.met"))
  sink("/dev/null");
  suppressMessages({
    mae <- createMAE(exp = gene.exp, met = dna.met, 
                     sampleMap = sampleMap, colData = sample.info, genome = "hg19") 
  })
  expect_equal(metadata(mae)$genome,"hg19")
  expect_false(metadata(mae)$TCGA)
  expect_equal(dim(getExp(mae)),c(2,2))
  expect_equal(dim(getMet(mae)),c(2,2))
  expect_equal(assay(getMet(mae)),as.matrix(dna.met))
  expect_equal(assay(getExp(mae)),as.matrix(gene.exp))
  expect_equal(colData(mae),sample.info)
  expect_true(all(sampleMap(mae)$assay %in% c("DNA methylation","Gene expression")))
  expect_true(all(c("external_gene_name","ensembl_gene_id") %in% colnames(values(getExp(mae)))))
})

test_that("The creation of a using Summarized Experiment objects and TCGA data", {
  # NON TCGA example: matrices has diffetrent column names
  gene.exp <- DataFrame(sample1.exp = c("ENSG00000141510"=2.3,"ENSG00000171862"=5.4),
                        sample2.exp = c("ENSG00000141510"=1.6,"ENSG00000171862"=2.3)
  )
  dna.met <- DataFrame(sample1.met = c("cg14324200"=0.5,"cg23867494"=0.1),
                       sample2.met =  c("cg14324200"=0.3,"cg23867494"=0.9)
  )
  sample.info <- DataFrame(sample.type = c("Normal", "Tumor"))
  rownames(sample.info) <- c("sample1","sample2")
  sampleMap <- DataFrame(primary = c("sample1","sample1","sample2","sample2"), 
                         colname = c("sample1.exp","sample1.met","sample2.exp","sample2.met"))
  sink("/dev/null");
  suppressMessages({
    mae <- createMAE(exp = gene.exp, met = dna.met, 
                     sampleMap = sampleMap, colData = sample.info, genome = "hg19") 
  })
  expect_equal(metadata(mae)$genome,"hg19")
  expect_false(metadata(mae)$TCGA)
  expect_equal(dim(getExp(mae)),c(2,2))
  expect_equal(dim(getMet(mae)),c(2,2))
  expect_equal(assay(getMet(mae)),as.matrix(dna.met))
  expect_equal(assay(getExp(mae)),as.matrix(gene.exp))
  expect_equal(colData(mae),sample.info)
  expect_true(all(sampleMap(mae)$assay %in% c("DNA methylation","Gene expression")))
  expect_true(all(c("external_gene_name","ensembl_gene_id") %in% colnames(values(getExp(mae)))))
  
})

test_that("The creation of a using Summarized Experiment objects and TCGA data", {
  
  # TCGA example using TCGAbiolinks
  # Testing creating MultyAssayExperiment object
  # Load library
  # Consisering it is TCGA and SE
  data(elmer.data.example, envir = environment())
  sink("/dev/null");
  suppressMessages({
    mae <- createMAE(exp = getExp(data), 
                     met =  getMet(data), 
                     TCGA = TRUE, genome = "hg19")
  })
  expect_equal(metadata(mae)$genome,"hg19")
  expect_true(metadata(mae)$TCGA)
  expect_true(all(sampleMap(mae)$assay %in% c("DNA methylation","Gene expression")))
  expect_true(all(c("external_gene_name","ensembl_gene_id") %in% colnames(values(getExp(mae)))))
  expect_equal(dim(getMet(mae)),c(101,234))
  expect_equal(dim(getExp(mae)),c(1026,234))
  
  sink("/dev/null");
  suppressMessages({
    mae <- createMAE(exp = getExp(data), 
                     met =  getMet(data), 
                     TCGA = TRUE, 
                     genome = "hg38")
  })
  expect_equal(metadata(mae)$genome,"hg38")
  expect_true(metadata(mae)$TCGA)
  expect_true(all(sampleMap(mae)$assay %in% c("DNA methylation","Gene expression")))
  expect_true(all(c("external_gene_name","ensembl_gene_id") %in% colnames(values(getExp(mae)))))
  expect_equal(dim(getMet(mae)),c(101,234))
  expect_equal(dim(getExp(mae)),c(1026,234))
  
  # Consisering it is TCGA and not SE
  sink("/dev/null");
  suppressMessages({
    mae <- createMAE(exp = assay(getExp(data)), met =  assay(getMet(data)), 
                     TCGA = TRUE, genome = "hg19")
  })
  expect_equal(metadata(mae)$genome,"hg19")
  sink("/dev/null");
  suppressMessages({
    mae <- createMAE(exp = assay(getExp(data)), 
                     met =  assay(getMet(data)), 
                     TCGA = TRUE, genome = "hg38")
  })
  expect_equal(metadata(mae)$genome,"hg38")
  
  # Consisering it is not TCGA and SE
  # DNA methylation and gene expression Objects should have same sample names in columns
  not.tcga.exp <- assay(getExp(data)) 
  colnames(not.tcga.exp) <- substr(colnames(not.tcga.exp),1,15)
  not.tcga.met <- assay(getMet(data)) 
  colnames(not.tcga.met) <- substr(colnames(not.tcga.met),1,15)
  
  phenotype.data <- data.frame(row.names = colnames(not.tcga.exp), 
                               samples = colnames(not.tcga.exp), 
                               group = c(rep("group1", length(colnames(not.tcga.exp)) / 2 ), 
                                         rep("group2", length(colnames(not.tcga.exp)) /2 )))
  sink("/dev/null");
  suppressMessages({
    mae <- createMAE(exp = not.tcga.exp, met =  not.tcga.met,
                     TCGA = FALSE, genome = "hg19", colData = phenotype.data)
  })
  
})


test_that("Number of probes in MAE matches the distal probes", {
  library(TCGAbiolinks)
  library(dplyr)
  gcimp.samples <- TCGAquery_subtype("lgg") %>% dplyr::filter(base::grepl("G-CIMP",Supervised.DNA.Methylation.Cluster,ignore.case = T))
  #-----------------------------------
  # 2 - Get data
  # ----------------------------------
  #-----------------------------------
  # 2.1 - DNA Methylation
  # ----------------------------------
  query <- GDCquery(project = "TCGA-LGG", 
                    data.category = "DNA Methylation",
                    platform = "Illumina Human Methylation 450", 
                    barcode =  gcimp.samples$patient[1:3])
  GDCdownload(query)
  met <- GDCprepare(query, save = FALSE)
  #-----------------------------------
  # 2 - Expression
  # ----------------------------------
  query <- GDCquery(project = "TCGA-LGG", 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - FPKM-UQ",
                    barcode =   gcimp.samples$patient[1:3])
  GDCdownload(query)
  exp <- GDCprepare(query, save = FALSE)
  for(genome in c("hg38","hg19")){
    distal.probe <- get.feature.probe(genome = genome, met.platform = "450K")
    mae <- createMAE(exp           = exp, 
                     met           = met,
                     genome        = genome, 
                     met.platform  = "450K", 
                     linearize.exp = TRUE, 
                     met.na.cut    = 1.1,
                     filter.probes = distal.probe, 
                     TCGA          = TRUE) 
    
    expect_equal(length(distal.probe),nrow(getMet(mae)))
    expect_equal(metadata(mae)$genome,genome)
    expect_true("TN" %in% colnames(colData(mae)))
  }
  unlink("GDCdata",recursive = TRUE, force = TRUE)
})
