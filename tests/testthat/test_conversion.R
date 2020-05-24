context("Test DESeq2")

test_that("Conversion to SummarizedExperiment works", {
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    obj <- readRDS("./obj.Rds")

    dds <- ThinDataToDESeqDataSet(obj = obj)
    dds2 <- DESeq2::DESeq(object = dds)
  }
})
