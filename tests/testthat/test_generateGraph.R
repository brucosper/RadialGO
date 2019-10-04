#test_countNonExon.R

context("radialgo")

# ==== BEGIN SETUP AND PREPARE =================================================
#

load("../../inst/data/go_ids.RData")
load("../../inst/data/testNodeDF.RData")
load("../../inst/data/testEdgeDF.RData")
gr <- generateGraph(go_ids, 5, 0.2)
ndf <- radialgo::get_node_df(gr)
edf <- radialgo::get_edge_df(gr)
#
# ==== END SETUP AND PREPARE ===================================================


test_that("Sample input produces the expected output",  {
  expect_equal(get_node_df(gr), ndf)
  expect_equal(get_edge_df(gr), edf)
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persitent construct that the test has created, except for
# stuff in tempdir().
#
rm(go_ids)
rm(ndf)
rm(edf)

# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
