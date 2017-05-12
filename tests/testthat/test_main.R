

context("Main and helper functions")

test_that("RepList function works", {
  param = MmapprParam(refGenome = new("GmapGenome"), "yo", "hi", VEPParam())
  expect_equal(length(RepList(param, 5)), 5)
  expect_equal(length(RepList(param, 0)), 0)
  expect_equal(length(RepList(param, 1)), 1)
  expect_equal(length(RepList(param, 2)), 2)
  expect_error(length(RepList(param, -1)))
  x <- RepList(param, 2)
  expect_identical(x[[1]], x[[2]])
})
