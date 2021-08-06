test_that("summarizeNeighborhood function works", {
   library(cytomapper)
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)

    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "CellType", method = "classic",
                                                   colPairName = "knn_interaction_graph"))     
})
