#test-generateGraph.R

context("radialgo")

# ==== BEGIN SETUP AND PREPARE =================================================
#

testNodeData <- system.file("testdata", "testNodeDF.RData", package="radialgo")
testEdgeData <- system.file("testdata", "testEdgeDF.RData", package="radialgo")
go_ids_loc <- system.file("testdata", "test.csv", package="radialgo")

go_ids <- read.csv(go_ids_loc)
load(testNodeData)
load(testEdgeData)
gr <- generateGraph(go_ids, 5, 0.2)
ndf <- DiagrammeR::get_node_df(gr)
edf <- DiagrammeR::get_edge_df(gr)
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
rm(ndf)
rm(edf)

# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
