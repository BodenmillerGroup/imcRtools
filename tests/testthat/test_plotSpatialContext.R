test_that("plotSpatialContext function works", {
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
  
  ## Plot spatial context - basic tests
  expect_silent(p <- plotSpatialContext(sce, group_by = "ImageNb"))
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_size_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_size_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_label_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_label_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_label_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  #return data - tests
  p <- plotSpatialContext(sce, group_by = "ImageNb", return_data = TRUE) 
  expect_type(p, "list")
  expect_equal(names(p), c("edges", "vertices"))
  #metadata entry and vertices comp
  sce_fil <- filterSpatialContext(sce, group_by = "ImageNb", group_threshold = 0)
  expect_equal(metadata(sce_fil)$filterSpatialContext, p$vertices[,1:3])
  
  #Errors
  expect_error(plotSpatialContext(colData(sce)),
               regexp = "'object' needs to be a SingleCellExperiment object.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, entry = "spatialcontext"),
               regexp = "'entry' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_color_by = "NAME"),
               regexp = "'node_color_by' has to be one off 'name','n_cells','n_group'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_size_by = "name"),
               regexp = "'node_size_by' has to be 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_color_fix = factor("black")),
               regexp = "'node_color_fix' has to be a character'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_size_fix = 3),
               regexp = "'node_size_fix' has to be a character'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_label_repel = "TRUE"),
               regexp = "'node_label_repel' has to be logical'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_label_color_by = "NAME"),
               regexp = "'node_label_color_by' has to be one off 'name','n_cells','n_group'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_label_color_fix = factor("black")),
               regexp = "'node_label_color_fix' has to be a character'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", draw_edges = "TRUE"),
               regexp = "'draw_edges' has to be logical'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", edge_color_fix = factor("black")),
               regexp = "'edge_color_fix' has to be a character'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", return_data = "TRUE"),
               regexp = "'return_data' has to be logical'.",
               fixed = TRUE)
}
)