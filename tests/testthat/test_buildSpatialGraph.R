test_that("buildSpatialGraph function works", {
    library(cytomapper)
    
    data("pancreasImages")
    data("pancreasMasks")

    cur_sce <- measureObjects(pancreasMasks, pancreasImages, img_id = "ImageNb")
    
    buildSpatialGraph(cur_sce, )
})
