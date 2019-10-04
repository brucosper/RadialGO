#test_countNonExon.R

context("radialgo")

# ==== BEGIN SETUP AND PREPARE =================================================
#

load("data/go_ids.RData")
load("data/testOutput")


#
# ==== END SETUP AND PREPARE ===================================================


test_that("Sample input prouces the expected output",  {
  expect_equal(generateGraph(go_ids, 5, 0.2), graph)
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persitent construct that the test has created, except for
# stuff in tempdir().
#
rm(go_ids)
rm(graph)

# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
