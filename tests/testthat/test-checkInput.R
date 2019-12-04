#test-checkInput.R

test_that("input check works for invalid data frame", {
  expect_equal(checkInput(data.frame(test1=c(1,2,3), test2=c(4,5,6)), 5, 0.5), FALSE)
  expect_warning(checkInput(data.frame(test1=c(1,2,3), test2=c(4,5,6)), 5, 0.5), "Invalid data frame. Please make sure to name the columns 'GO ID' and 'P value'.")

})
test_that("input check works for invalid top value", {
  expect_equal(checkInput(data.frame(GO.ID="GO:0008150", P.value=0.005), "a", 0.5), FALSE)
  expect_warning(checkInput(data.frame(GO.ID="GO:0008150", P.value=0.005), "a", 0.5), "Invalid top value.")
})
test_that("input check works for invalid cutoff value", {
  expect_equal(checkInput(data.frame(GO.ID="GO:0008150", P.value=0.005), 5, "a"), FALSE)
  expect_warning(checkInput(data.frame(GO.ID="GO:0008150", P.value=0.005), 5, "a"), "Invalid cutoff value.")

})

# [END]
