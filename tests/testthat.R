library(testthat)
library(MMAPPR2)

BiocParallel::register(BiocParallel::SerialParam())

testthat::test_check("MMAPPR2")
