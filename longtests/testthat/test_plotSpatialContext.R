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
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", draw_edges = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(sce, group_by = "ImageNb", node_color_by = "name",
                          node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce$spatial_context)),"_"),length))
  
  # filtered entries - tests
  sce_fil <- filterSpatialContext(sce, group_by = "ImageNb", group_threshold = 2)
  
  p <- plotSpatialContext(sce_fil, group_by = "ImageNb", entry = "spatial_context_filtered")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(sce_fil$spatial_context_filtered)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(sce_fil$spatial_context_filtered)),"_"),length))
  
  # return data - tests
  p <- plotSpatialContext(sce, group_by = "ImageNb", return_data = TRUE) 
  expect_type(p, "list")
  expect_equal(names(p), c("edges", "vertices"))
  expect_equal(p$edges[,1], c("1", "1", "1_2", "1_3", "2", "2", "2_3", "3", "3"))
  expect_equal(p$edges[,2], c("1_2", "1_3", "1_2_3", "1_2_3", "1_2", "2_3", "1_2_3", 
                              "1_3", "2_3"))
  expect_equal(p$vertices[,1], c("1", "1_2", "1_2_3", "1_3", "2", "2_3", "3"))
  expect_equal(p$vertices[,2], c(87L, 71L, 16L, 90L, 55L, 29L, 14L))
  
  #metadata entry and vertices comp
  sce_fil <- filterSpatialContext(sce, group_by = "ImageNb", group_threshold = 0)
  expect_equal(metadata(sce_fil)$filterSpatialContext, p$vertices[,1:3])
  
  # aggregatedNeighbors colnames as characters
  cur_sce <- sce
  colnames(cur_sce$aggregatedNeighborhood) <- c("B_CN","T_CN","DC_CN")
  expect_silent(cur_sce <- detectSpatialContext(cur_sce, entry = "aggregatedNeighborhood",
                                                  threshold = 0.9))
  expect_silent(cur_sce <- filterSpatialContext(cur_sce, group_by = "ImageNb", 
                                                  group_threshold = 2))
  
  p <- plotSpatialContext(cur_sce, group_by = "ImageNb")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(cur_sce$spatial_context)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(cur_sce$spatial_context)),"_"),length))
  
  p <- plotSpatialContext(cur_sce, group_by = "ImageNb", entry = "spatial_context_filtered")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(cur_sce$spatial_context_filtered)))
  expect_equal(p$data$length, sapply(str_split(sort(unique(cur_sce$spatial_context_filtered)),"_"),length))
  
  #Errors
  expect_error(plotSpatialContext(colData(sce)),
               regexp = "'object' needs to be a SingleCellExperiment object.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, entry = "spatialcontext"),
               regexp = "'entry' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "Image"),
               regexp = "'group_by' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_color_by = "NAME"),
               regexp = "'node_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
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
               regexp = "'node_label_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
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
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_label_repel = FALSE, 
                                  node_label_color_by = "name"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined when node_label_repel == FALSE",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_label_repel = FALSE, 
                                  node_label_color_by = "name"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined when node_label_repel == FALSE",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_color_by = "name", 
                                  node_color_fix = "blue"),
               regexp = "'node_color_by' and 'node_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_label_color_by = "name", 
                                  node_label_color_fix = "blue"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_color_by = "name", 
                                  node_label_color_by = "n_cells"),
               regexp = "'node_label_color_by' and 'node_color_by' have to be identical.",
               fixed = TRUE)
  
  expect_error(plotSpatialContext(sce, group_by = "ImageNb", node_size_by = "n_group", 
                                  node_size_fix = "22"),
               regexp = "'node_size_by' and 'node_size_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  
  #spatial context is NA vector
  sce$spatial_context <- NA
  expect_error(plotSpatialContext(sce, group_by = "ImageNb"), 
               regexp = "the data frame should contain at least two columns", 
               fixed = TRUE)
  
}
)