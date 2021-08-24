test_that("plotSpatial function works", {
    library(cytomapper)
    library(ggraph)
    data("pancreasSCE")
    
    cur_sce <- pancreasSCE

    # SingleCellExperiment
    expect_silent(p <- plotSpatial(cur_sce, img_id = "ImageNb"))
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
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "PIN", assay_type = "counts")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$PIN, counts(pancreasSCE)["PIN",], check.attributes = FALSE)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "PIN", assay_type = "exprs")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$PIN, assay(pancreasSCE, "exprs")["PIN",], check.attributes = FALSE)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "CDH", assay_type = "exprs")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$CDH, assay(pancreasSCE, "exprs")["CDH",], check.attributes = FALSE)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_shape_by = "CellType")
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    expect_equal(p$data$x, pancreasSCE$Pos_X)
    expect_equal(p$data$y, pancreasSCE$Pos_Y)
    expect_equal(p$data$ImageNb, pancreasSCE$ImageNb)
    expect_equal(p$data$CellType, pancreasSCE$CellType)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", node_shape_by = "Area")
    expect_s3_class(p, "ggraph")
    expect_warning(print(p))
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
                     edge_color_by = "CellType", node_color_by = "PIN",
                     assay_type = "exprs", node_size_fix = 5)
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
                     edge_color_by = "Pattern", node_color_by = "Pattern",
                     nodes_first = FALSE)
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
    
    p2 <- plotSpatial(cur_sce, img_id = "ImageNb", 
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     directed = FALSE)
    expect_s3_class(p2, "ggraph")
    expect_silent(print(p2))
    expect_equal(p2$data$x, pancreasSCE$Pos_X)
    expect_equal(p2$data$y, pancreasSCE$Pos_Y)
    expect_equal(p2$data$ImageNb, pancreasSCE$ImageNb)
    
    cur_graph <- igraph::as.igraph(attributes(p$data)$graph)
    cur_graph <- as.undirected(cur_graph)
    cur_graph_2 <- igraph::as.igraph(attributes(p2$data)$graph)
    
    cur_edges <- as_edgelist(cur_graph)
    cur_edges_2 <- as_edgelist(cur_graph_2)

    cur_edges <- cur_edges[order(paste(cur_edges[,1], cur_edges[,2])),]
    cur_edges_2 <- cur_edges_2[order(paste(cur_edges_2[,1], cur_edges_2[,2])),]
    
    expect_equal(cur_edges, cur_edges_2)
    
    # Delaunay
    cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "delaunay")
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "delaunay_interaction_graph",
                     edge_width_fix = 0.1)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p2 <- plotSpatial(cur_sce, img_id = "ImageNb", 
                      draw_edges = TRUE, colPairName = "delaunay_interaction_graph",
                      directed = FALSE)
    expect_s3_class(p2, "ggraph")
    expect_silent(print(p2))
    expect_equal(p2$data$x, pancreasSCE$Pos_X)
    expect_equal(p2$data$y, pancreasSCE$Pos_Y)
    expect_equal(p2$data$ImageNb, pancreasSCE$ImageNb)
    
    cur_graph <- igraph::as.igraph(attributes(p$data)$graph)
    cur_graph <- as.undirected(cur_graph)
    cur_graph_2 <- igraph::as.igraph(attributes(p2$data)$graph)
    
    cur_edges <- as_edgelist(cur_graph)
    cur_edges_2 <- as_edgelist(cur_graph_2)
    
    cur_edges <- cur_edges[order(paste(cur_edges[,1], cur_edges[,2])),]
    cur_edges_2 <- cur_edges_2[order(paste(cur_edges_2[,1], cur_edges_2[,2])),]
    
    expect_equal(cur_edges, cur_edges_2)
    
    # Expansion
    cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "expansion", threshold = 15)
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", draw_edges = TRUE, colPairName = "expansion_interaction_graph",
                     edge_width_fix = 0.1)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p2 <- plotSpatial(cur_sce, img_id = "ImageNb", 
                      draw_edges = TRUE, colPairName = "expansion_interaction_graph",
                      directed = FALSE)
    expect_s3_class(p2, "ggraph")
    expect_silent(print(p2))
    expect_equal(p2$data$x, pancreasSCE$Pos_X)
    expect_equal(p2$data$y, pancreasSCE$Pos_Y)
    expect_equal(p2$data$ImageNb, pancreasSCE$ImageNb)
    
    cur_graph <- igraph::as.igraph(attributes(p$data)$graph)
    cur_graph <- as.undirected(cur_graph)
    cur_graph_2 <- igraph::as.igraph(attributes(p2$data)$graph)
    
    cur_edges <- as_edgelist(cur_graph)
    cur_edges_2 <- as_edgelist(cur_graph_2)
    
    cur_edges <- cur_edges[order(paste(cur_edges[,1], cur_edges[,2])),]
    cur_edges_2 <- cur_edges_2[order(paste(cur_edges_2[,1], cur_edges_2[,2])),]
    
    expect_equal(cur_edges, cur_edges_2)
    
    
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
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", node_color_by = "CellType",
                     flip_y = FALSE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", node_color_by = "CellType",
                     flip_x = TRUE)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageName", node_color_by = "CellType",
                     flip_x = TRUE, flip_y = FALSE)
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
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE, 
                     colPairName = "knn_interaction_graph",
                     arrow = arrow(angle = 10, length = unit(0.05, "inch"), type = "closed"),
                     end_cap = circle(0.3, "cm"),
                     directed = TRUE, scales = "free", nrow = 1)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE, 
                     colPairName = "knn_interaction_graph",
                     arrow = arrow(angle = 10, length = unit(0.05, "inch"), type = "closed"),
                     end_cap = circle(0.3, "cm"),
                     directed = TRUE, scales = "free", ncol = 1)
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))

    ## Modify plots
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE, 
                     colPairName = "knn_interaction_graph", node_color_by = "CellType",
                     node_shape_by = "ImageName", node_size_by = "Area", edge_color_by = "Pattern",
                     edge_width_by = "Area") + 
        scale_edge_color_brewer(palette = "Set1") +
        scale_color_manual(values = c(celltype_A = "yellow",
                                      celltype_C = "green",
                                      celltype_B = "dark blue")) +
        scale_shape_manual(values = c(E34_imc.tiff = 10,
                                      G01_imc.tiff = 12,
                                      J02_imc.tiff = 14)) +
        scale_edge_width(range = c(0.1,2)) +
        scale_size(range = c(4, 12))
    expect_s3_class(p, "ggraph")
    expect_silent(print(p))
    
    # SpatialExperiment
    library(SpatialExperiment)
    cur_spe <- SpatialExperiment(assays = list(counts = counts(cur_sce)),
                                 sample_id = cur_sce$ImageName) 
    colData(cur_spe) <- colData(cur_sce)
    colPairs(cur_spe) <- colPairs(cur_sce)
    
    # Error
    expect_error(p <- plotSpatial(cur_spe, img_id = "ImageNb", draw_edges = TRUE, 
                     colPairName = "knn_interaction_graph", node_color_by = "CellType",
                     node_shape_by = "ImageName", node_size_by = "Area", edge_color_by = "Pattern",
                     edge_width_by = "Area"), regex = "'coords' not in spatialCoords(object).",
                 fixed = TRUE)
    
    spatialCoords(cur_spe) <- as.matrix(colData(cur_sce)[,c("Pos_X", "Pos_Y")])
    p1 <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE, 
                      colPairName = "knn_interaction_graph", node_color_by = "CellType",
                      node_shape_by = "ImageName", node_size_by = "Area", edge_color_by = "Pattern",
                      edge_width_by = "Area")
    p2 <- plotSpatial(cur_spe, img_id = "ImageNb", draw_edges = TRUE, 
                     colPairName = "knn_interaction_graph", node_color_by = "CellType",
                     node_shape_by = "ImageName", node_size_by = "Area", edge_color_by = "Pattern",
                     edge_width_by = "Area")
    expect_s3_class(p2, "ggraph")
    expect_silent(print(p2))
    expect_equal(p$data, p2$data, check.attributes = FALSE)
    expect_equal(p$layers, p2$layers, check.attributes = FALSE)
    expect_equal(p$mapping, p2$mapping, check.attributes = FALSE)
    expect_equal(p$theme, p2$theme, check.attributes = FALSE)
    expect_equal(p$coordinates, p2$coordinates, check.attributes = FALSE)
    expect_equal(p$plot_env, p2$plot_env, check.attributes = FALSE)
    expect_equal(p$labels, p2$labels, check.attributes = FALSE)
    
    # Fail
    expect_error(plotSpatial("test"), "'object' not of type 'SingleCellExperiment'.")
    expect_error(plotSpatial(cur_sce), 'argument "img_id" is missing, with no default')
    expect_error(plotSpatial(cur_sce, img_id = "test"), "'img_id' not in colData(object).",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = 1), "'img_id' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", coords = 1), 
                 "'coords' must be a character vector of length 2.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", coords = c(1, 2)), 
                 "'coords' must be a character vector of length 2.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", coords = c("Pos_X", "test")), 
                 "'coords' not in colData(object).",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_spe, img_id = "ImageNb", coords = c("Pos_X", "test")), 
                 "'coords' not in spatialCoords(object).",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = 1), 
                 "'node_color_by' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = c("test", "test2")), 
                 "'node_color_by' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "test"), 
                 "'node_color_by' not in colData(object) or rownames(object).",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "PIN"), 
                 "When coloring nodes by marker expression, please specify 'assay_type'.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "PIN", assay_type = 1), 
                 "'assay_type' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "PIN", assay_type = "test"), 
                 "'assay_type' not an assay in object.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_shape_by = c("test", "test2")), 
                 "'node_shape_by' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_shape_by = "test"), 
                 "'node_shape_by' not in colData(object).",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_size_by = c("test", "test2")), 
                 "'node_size_by' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", node_size_by = "test"), 
                 "'node_size_by' not in colData(object).",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = "test"), 
                 "'draw_edges' must be a single logical",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = c("test", "test2")), 
                 "'draw_edges' must be a single logical",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE), 
                 "Please specify the name of the column pairing via 'colPairName'.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = "test"), 
                 "No column pairing with name test found.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = 1), 
                 "'colPairName' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = c("test", "test2")), 
                 "'colPairName' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = "knn_interaction_graph",
                             edge_color_by = "test2"), 
                 "'edge_color_by' not in 'colData(object)' or in 'mcols(colPair(object, colPairName))'.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = "knn_interaction_graph",
                             edge_color_by = c("test", "test2")), 
                 "'edge_color_by' must be a single string.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = "knn_interaction_graph",
                             edge_width_by = "test2"), 
                 "'edge_width_by' not in 'colData(object)' or in 'mcols(colPair(object, colPairName))'.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = "knn_interaction_graph",
                             edge_width_by = c("test", "test2")), 
                 "'edge_width_by' must be a single string.",
                 fixed = TRUE)
    
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = "knn_interaction_graph",
                             directed = "test"), 
                 "'directed' must be a single logical",
                 fixed = TRUE)
    
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = "knn_interaction_graph",
                             arrow = "test"), 
                 "'arrow' must be of class grid::arrow.",
                 fixed = TRUE)
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                             colPairName = "knn_interaction_graph",
                             arrow = grid::arrow(),
                             end_cap = "test"), 
                 "'end_cap' must be of type 'geometry'.",
                 fixed = TRUE)

    expect_error(plotSpatial(cur_sce, img_id = "ImageNb",
                             nodes_first = "test"), 
                 "'nodes_first' must be a single logical",
                 fixed = TRUE)
        
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb",
                             ncols = "test"), 
                 "'ncols' must be a single numeric",
                 fixed = TRUE)
    
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb",
                             nrows = "test"), 
                 "'nrows' must be a single numeric",
                 fixed = TRUE)
    
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb",
                             scales = "test"), 
                 "'scales' should be one of 'fixed', 'free_x', 'free_y', 'free'.",
                 fixed = TRUE)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb",
                     node_color_fix = "test")
    
    expect_error(print(p), 
                 "Unknown colour name: test",
                 fixed = TRUE)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb",
                     node_shape_fix = "test")
    
    expect_error(print(p), 
                 "Can't find shape name:\n* 'test'",
                 fixed = TRUE)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb",
                     node_size_fix = "test")
    
    expect_error(print(p), 
                 "non-numeric argument to binary operator",
                 fixed = TRUE)
    
    p <- plotSpatial(cur_sce, img_id = "ImageNb",
                     draw_edges = TRUE, colPairName = "knn_interaction_graph",
                     edge_color_fix = "test")
    
    expect_error(print(p), 
                 "Unknown colour name: test",
                 fixed = TRUE)
    
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb",
                             flip_x = "test"), 
                 "'flip_x' must be a single logical",
                 fixed = TRUE)
    
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb",
                             flip_x = c(1, 2)), 
                 "'flip_x' must be a single logical",
                 fixed = TRUE)
    
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb",
                             flip_y = "test"), 
                 "'flip_y' must be a single logical",
                 fixed = TRUE)
    
    expect_error(plotSpatial(cur_sce, img_id = "ImageNb",
                             flip_y = c(1, 2)), 
                 "'flip_y' must be a single logical",
                 fixed = TRUE)
    
})
