test_that("filterSpatialContext function works", {
  set.seed(22)
  library(cytomapper)
  data(pancreasSCE)
  
  ## 1. Cellular neighborhood (CN)
  sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb",
                           type = "knn",
                           name = "knn_cn_graph",
                           k = 5)
  
  sce <- aggregateNeighbors(sce, colPairName = "knn_cn_graph",
                            aggregate_by = "metadata",
                            count_by = "CellType",
                            name = "aggregatedCellTypes")
  
  cur_cluster <- kmeans(sce$aggregatedCellTypes, centers = 3)
  sce$cellular_neighborhood <- factor(cur_cluster$cluster)
  
  ## 2. Spatial context (SC)
  sce <- buildSpatialGraph(sce, img_id = "ImageNb",
                           type = "knn",
                           name = "knn_sc_graph",
                           k = 15)
  
  sce <- aggregateNeighbors(sce, colPairName = "knn_sc_graph",
                            aggregate_by = "metadata",
                            count_by = "cellular_neighborhood",
                            name = "aggregatedNeighborhood")
  
  # Detect spatial context
  sce <- detectSpatialContext(sce, entry = "aggregatedNeighborhood",
                              threshold = 0.9)
  
  ## Filter spatial context - tests
  # basics
  expect_silent(cur_sce <- filterSpatialContext(sce, group_by = "ImageNb", group_threshold = 2))
  
  expect_equal(names(colData(cur_sce)), 
               c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb", 
                 "CellNb", "MaskName", "Pattern", "aggregatedCellTypes", "cellular_neighborhood", 
                 "aggregatedNeighborhood", "spatial_context", "spatial_context_filtered"))
  
  expect_s4_class(cur_sce , class = "SingleCellExperiment") 
  expect_type(cur_sce$spatial_context_filtered, "character")
  expect_equal(length(cur_sce$spatial_context), length(cur_sce$spatial_context_filtered))
  
  # check that detected spatial_context_filtered remain the same
  expect_equal(cur_sce$spatial_context_filtered[200:210], c("3", NA, NA, "1_3", "1_3", "1_3", "1_3", NA, NA, "3", "1_3"))
  
  # change filtering
  expect_silent(cur_sce2 <- filterSpatialContext(sce, group_by = "ImageNb", cells_threshold = 15))
  
  expect_false(identical(cur_sce$spatial_context_filtered, cur_sce2$spatial_context_filtered))
  expect_equal(cur_sce$spatial_context, cur_sce2$spatial_context)
  
  # filtering with 0
  expect_silent(cur_sce3 <- filterSpatialContext(sce, group_by = "ImageNb", 
                                                 group_threshold = 0, 
                                                 cells_threshold = 0))
  expect_equal(cur_sce3$spatial_context, cur_sce3$spatial_context_filtered)
  
  # filtering with values higher than input
  expect_silent(cur_sce4 <- filterSpatialContext(sce, group_by = "ImageNb", 
                                                 group_threshold = 4, 
                                                 cells_threshold = 100))
  expect_equal(as.character(cur_sce4$spatial_context_filtered), 
               as.character(rep(NA, length(cur_sce4$spatial_context_filtered))))
  
  # manual filtering should give the same results
  cur_anno <- colData(sce) %>% as.data.frame() %>% select("spatial_context", "ImageNb") %>% 
    table() %>% as.data.frame()
  cur_anno <- data.frame(spatial_context = unique(cur_anno$spatial_context) %>% unfactor(), 
                     n_cells = cur_anno %>% group_by_at("spatial_context") %>% 
                       summarise(sum = sum(Freq)) %>% pull(sum), 
                     n_group = cur_anno %>% group_by_at("spatial_context") %>% 
                       filter(Freq != 0) %>% count() %>% pull(n)
  )
  
  manual_selected <- cur_anno %>% filter(n_group >= 2 & n_cells >= 15) %>% 
    pull(spatial_context) %>% sort()
  
  expect_silent(cur_sce5 <- filterSpatialContext(sce, group_by = "ImageNb", 
                                                 group_threshold = 2, 
                                                 cells_threshold = 15))
  
  function_selected <- unique(cur_sce5$spatial_context_filtered) %>% sort()

  expect_equal(manual_selected, function_selected)
  
  #metadata entry check
  manual_anno <- cur_anno[cur_anno$n_group >= 2,]
  expect_equal(metadata(cur_sce)$filterSpatialContext, manual_anno)
  
  #Errors
  expect_error(filterSpatialContext(colData(sce)),
               regexp = "'object' needs to be a SingleCellExperiment object.",
               fixed = TRUE)

  expect_error(filterSpatialContext(sce, entry = "spatialcontext"),
               regexp = "'entry' not in 'colData(object)'.",
               fixed = TRUE)

  expect_error(filterSpatialContext(sce, group_by = "ImageNb", cells_threshold = "10"),
               regexp = "'cells_threshold' needs to be a single numeric.",
               fixed = TRUE)
  
  expect_error(filterSpatialContext(sce, group_by = "Image", group_threshold = "2"),
               regexp = "'group_by' not in 'colData(object)'.",
               fixed = TRUE)

  expect_error(filterSpatialContext(sce, group_by = "ImageNb", group_threshold = "2"),
               regexp = "'group_threshold' needs to be a single numeric.",
               fixed = TRUE)
  
  expect_error(filterSpatialContext(sce, group_by = "ImageNb"),
               regexp = "One of 'group_threshold' and 'cells_threshold' has to be defined.",
               fixed = TRUE)
  
  expect_error(filterSpatialContext(sce, group_by = "ImageNb", cells_threshold = 10, 
                                    name = c("spatial","context","filtered")),
               regexp = "'name' has to be a single character'.",
               fixed = TRUE)
}
)
  