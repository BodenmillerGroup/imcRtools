test_that("buildSpatialGraph function works", {
    library(cytomapper)
    data("pancreasSCE")
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "expansion", threshold = 10)
})
