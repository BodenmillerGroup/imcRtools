test_that("detectSpatialContext function works", {
  library(cytomapper)
  data(pancreasSCE)
  
  ## 1. Cellular neighborhood (CN)
  sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                           type = "knn", 
                           name = "knn_cn_graph", 
                           k = 5)
  
  sce <- aggregateNeighbors(sce, colPairName = "knn_cn_graph", 
                            aggregate_by = "metadata", 
                            count_by = "CellType")
  
  cur_cluster <- kmeans(sce$aggregatedNeighbors, centers = 3)
  sce$cellular_neighborhood <- factor(cur_cluster$cluster)
  
  ## 2. Spatial context (SC)
  sce <- buildSpatialGraph(sce, img_id = "ImageNb", 
                           type = "knn", 
                           name = "knn_sc_graph", 
                           k = 15)
  
  sce <- aggregateNeighbors(sce, colPairName = "knn_sc_graph", 
                            aggregate_by = "metadata", 
                            count_by = "cellular_neighborhood")
  
  # Detect spatial context - tests
  expect_silent(cur_sce <- detectSpatialContext(sce, threshold = 0.9))
  expect_equal(names(colData(cur_sce)), 
               c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                 "CellNb", "MaskName", "Pattern", "aggregatedNeighbors", 
                 "cellular_neighborhood", "spatial_context"))
  expect_s4_class(cur_sce , class = "SingleCellExperiment") 
  expect_equal(length(cur_sce$spatial_context), length(cur_sce$cellular_neighborhood))
  
  # change threshold
  expect_silent(cur_sce_2 <- detectSpatialContext(sce, threshold = 0.5))
  expect_false(identical(cur_sce$spatial_context, cur_sce_2$spatial_context))
  
  #aggregatedNeighbors DataFrame contains a row with 0s
  cur_sce_2$aggregatedNeighbors[1,] <- 0
  expect_silent(cur_sce_3 <- detectSpatialContext(cur_sce_2, threshold = 0.9))
  expect_true(is.na(cur_sce_3$spatial_context[1]))
  
}
)