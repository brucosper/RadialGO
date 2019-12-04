#test-rebuildGOID.R

context("radialgo")

test_that("Sample input prouces the expected output",  {
  expect_equal(rebuildGOID(8150), "GO:0008150")
})

# [END]
