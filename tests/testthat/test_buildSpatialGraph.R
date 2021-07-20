test_that("buildSpatialGraph function works", {
    library(cytomapper)
    data("pancreasSCE")
    
    # Delauney
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "delaunay")
    expect_equal(colPairNames(cur_sce), "delaunay_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 1041)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                colPairName = "delaunay_interaction_graph")
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6))
    expect_equal(to(colPair(cur_sce))[10:20], c(3,  6,  7, 10,  2, 13, 24,  2,  5,  9, 12))
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "delaunay", directed = FALSE)
    expect_equal(length(colPair(cur_sce)), 2082)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                              colPairName = "delaunay_interaction_graph")
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4))
    expect_equal(to(colPair(cur_sce))[10:20], c(6, 11, 13, 50,  4,  6, 10, 12,  1,  2,  3))
    
    # KNN
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 5)
    expect_equal(colPairNames(cur_sce), "knn_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 1810)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                              colPairName = "knn_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                              colPairName = "knn_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4))
    expect_equal(to(colPair(cur_sce))[10:20], c(22,  4,  6, 10, 12, 17,  3,  7,  8, 10, 16))
    
    cur_ks <- table(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_ks <- rowSums(cur_ks)
    expect_true(all(cur_ks == 5))
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){which(x == 1)})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){which(x == 1)})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){which(x == 1)})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 5, directed = FALSE, mutual_edges = FALSE)
    expect_equal(colPairNames(cur_sce), "knn_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 1034)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5))
    expect_equal(to(colPair(cur_sce))[10:20], c(22,  4,  6, 10, 12, 17,  7,  8, 10, 16,  9))
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 5, directed = FALSE)
    expect_equal(colPairNames(cur_sce), "knn_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 2068)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4))
    expect_equal(to(colPair(cur_sce))[10:20], c(22,  4,  6, 10, 12, 17,  3,  7,  8, 10, 16))
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.directed(as.undirected(cur_ind_real))
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.directed(as.undirected(cur_ind_real))
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.directed(as.undirected(cur_ind_real))
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    # Other ks
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 10)
    expect_equal(colPairNames(cur_sce), "knn_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 3620)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    
    cur_ks <- table(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_ks <- rowSums(cur_ks)
    expect_true(all(cur_ks == 10))
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){which(x == 1)})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){which(x == 1)})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){which(x == 1)})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 10, directed = FALSE, mutual_edges = FALSE)
    expect_equal(colPairNames(cur_sce), "knn_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 2062)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 10, directed = FALSE)
    expect_equal(colPairNames(cur_sce), "knn_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 4124)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
    expect_equal(to(colPair(cur_sce))[10:20], c(29,  5,  9, 11, 13, 14, 19, 22, 24, 28, 31))
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.directed(as.undirected(cur_ind_real))
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.directed(as.undirected(cur_ind_real))
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.directed(as.undirected(cur_ind_real))
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    # Other algorithm
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 10)
    cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "knn", k = 10, name = "vptree",
                                 BNPARAM = VptreeParam())
    cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "knn", k = 10, name = "annoy",
                                 BNPARAM = AnnoyParam())
    cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "knn", k = 10, name = "hnsw",
                                 BNPARAM = HnswParam())
    
    expect_equal(colPairNames(cur_sce), 
                 c("knn_interaction_graph", "vptree", "annoy", "hnsw"))
    expect_equal(colPair(cur_sce, "knn_interaction_graph"),
                 colPair(cur_sce, "vptree"))
    
    # Expansion
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "expansion", threshold = 15)
    expect_equal(colPairNames(cur_sce), "expansion_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 2042)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "expansion_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "expansion_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6))
    expect_equal(to(colPair(cur_sce))[10:20], c(3, 7,  8, 10,  2, 11, 13, 19,  3,  9, 12))
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "expansion", threshold = 15, directed = FALSE)
    expect_equal(colPairNames(cur_sce), "expansion_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 2042)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "expansion_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "expansion_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6))
    expect_equal(to(colPair(cur_sce))[10:20], c(3, 7,  8, 10,  2, 11, 13, 19,  3,  9, 12))
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "expansion", threshold = 15, 
                                 directed = FALSE, mutual_edges = FALSE)
    expect_equal(colPairNames(cur_sce), "expansion_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 1021)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "expansion_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "expansion_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7))
    expect_equal(to(colPair(cur_sce))[10:20], c(7,  8, 10, 11, 13, 19,  9, 12,  8, 10, 16))
    
    # Parallelisation
    
    # SpatialExperiment
    
    # Fail
    
    
})
