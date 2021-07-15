test_that("plotSpatial function works", {
    library(cytomapper)
    data("pancreasSCE")
    
    cur_sce <- pancreasSCE

    # SingleCellExperiment
    p <- plotSpatial(cur_sce, img_id = "ImageNb")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$CellType, pancreasSCE$CellType)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "Area")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$Area, pancreasSCE$Area)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "Pattern")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$Pattern, pancreasSCE$Pattern)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_shape_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$CellType, pancreasSCE$CellType)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_shape_by = "Area")
    expect_s3_class(p, "ggraph")
    expect_warning(print(p), regexp = "The shape palette can deal with a maximum of 6 discrete values because more than 6 becomes difficult to discriminate; you
have 121. Consider specifying shapes manually if you must have them.")
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$Area, as.character(pancreasSCE$Area))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_shape_by = "Pattern")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$Pattern, as.character(pancreasSCE$Pattern))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_size_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_warning(print(p), regexp = "Using size for a discrete variable is not advised.")
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$CellType, pancreasSCE$CellType)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_size_by = "Area")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$Area, pancreasSCE$Area)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_size_by = "Pattern")
    expect_s3_class(p, "ggraph")
    expect_warning(print(p), regexp = "Using size for a discrete variable is not advised.")
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$Pattern, pancreasSCE$Pattern)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "CellType",
                     node_shape_by = "Pattern", node_size_by = "Area")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$Pattern, as.character(pancreasSCE$Pattern))
    expect_equal(p$data$CellType, pancreasSCE$CellType)
    expect_equal(p$data$Area, pancreasSCE$Area)
    
    cur_sce <- pancreasSCE
    cur_sce$ImageNb <- as.factor(cur_sce$ImageNb) 
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    cur_sce <- pancreasSCE
    cur_sce$ImageNb <- as.numeric(cur_sce$ImageNb) 
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", node_color_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", node_color_by = "CellType",
                     node_color_fix = "red")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName",
                     node_color_fix = "red")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName",
                     node_color_fix = 1)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", node_shape_by = "ImageNb",
                     node_shape_fix = "+")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName",
                     node_shape_fix = "+")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", node_shape_by = "ImageNb",
                     node_shape_fix = 21)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName",
                     node_shape_fix = 21)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", node_size_by = "Area",
                     node_size_fix = 5)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName",
                     node_size_fix = 5)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 3)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "CellType", node_color_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "Pattern", node_color_by = "Pattern")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "Area", node_color_by = "Area")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "Pattern")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "Area")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    mcols(colPair(cur_sce, "knn_interaction_graph"))$test <- runif(length(colPair(cur_sce, "knn_interaction_graph")))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "test")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    mcols(colPair(cur_sce, "knn_interaction_graph"))$test_2 <- letters[sample(x = 1:10,
                                            size = length(colPair(cur_sce, "knn_interaction_graph")),
                                            replace = TRUE)]
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "test_2")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "CellType",
                     edge_color_fix = "red")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_fix = "red")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_fix = 0.1, edge_width_by = "Area")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_fix = 0.1)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "CellType",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "CellType", node_color_by = "CellType",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "Pattern", node_color_by = "Pattern",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "Area", node_color_by = "Area",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "CellType",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "Pattern",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "Area",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "test",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))

    p <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "test_2",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_by = "CellType",
                     edge_color_fix = "red",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_fix = "red",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_fix = 0.1, edge_width_by = "Area",
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_fix = 0.1,
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_width_by = "CellType", edge_color_by = "Area", node_color_by = "CellType",
                     node_shape_by = "ImageNb", node_size_by = "Area",
                     directed = TRUE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     arrow = arrow(),
                     directed = TRUE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     arrow = arrow(angle = 10, length = unit(0.1, "inch"), type = "closed"),
                     directed = TRUE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     arrow = arrow(angle = 10, length = unit(0.1, "inch"), type = "closed"),
                     end_cap = circle(0.3, "cm"),
                     directed = TRUE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     arrow = arrow(angle = 10, length = unit(0.1, "inch"), type = "closed"),
                     end_cap = circle(0.3, "cm"),
                     directed = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    ## Subsetting
    cur_sce2 <- cur_sce[,cur_sce$Pattern]
    
    p <- plotSpatial(cur_sce2, img_id = "ImageNb", draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     arrow = arrow(angle = 10, length = unit(0.05, "inch"), type = "closed"),
                     end_cap = circle(0.3, "cm"),
                     directed = TRUE, scales = "fixed", node_size_fix = 0.1, edge_width_fix = 0.5)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce2, img_id = "ImageNb", draw_edges = TRUE, 
                     colPairName = "knn_interaction_graph", edge_width_by = "CellType", 
                     edge_color_by = "Area", node_color_by = "CellType",
                     node_shape_by = "ImageNb", node_size_by = "Area",)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce2[,cur_sce2$ImageNb == 2], img_id = "ImageNb", draw_edges = TRUE, 
                     colPairName = "knn_interaction_graph",
                     arrow = arrow(angle = 10, length = unit(0.05, "inch"), type = "closed"),
                     end_cap = circle(0.3, "cm"),
                     directed = TRUE, scales = "fixed", node_size_fix = 0.1, 
                     edge_width_fix = 0.5)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    

    ## Modify plots
    
    
    # SpatialExperiment
    
    # Fail
})
