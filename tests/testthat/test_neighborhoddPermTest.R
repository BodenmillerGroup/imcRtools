test_that("neighborhoodPermTest function works", {
    library(cytomapper)
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    ################################ classic ###################################
    set.seed(123)
    expect_silent(cur_out <- neighborhoodPermTest(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100))
    
    set.seed(123)
    expect_silent(cur_out_2 <- neighborhoodPermTest(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100))
    
    expect_identical(cur_out, cur_out_2)
    
    p <- BiocParallel::SnowParam(RNGseed = 123)
    expect_silent(cur_out_3 <- neighborhoodPermTest(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "classic",
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100,
                                                    BBPARAM = p))

    expect_silent(cur_out_3 <- neighborhoodPermTest(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "classic",
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100,
                                                    BBPARAM = p))
    
    expect_equal(cur_out_3, cur_out_4)
    
    

})
