test_that("buildSpatialGraph function works", {
    library(cytomapper)
    data("pancreasSCE")
    
    # Delauney
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "delaunay"))
    expect_equal(colPairNames(cur_sce), "delaunay_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 2082)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                colPairName = "delaunay_interaction_graph")
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4))
    expect_equal(to(colPair(cur_sce))[10:20], c(6, 11, 13, 50,  4,  6, 10, 12,  1,  2,  3))
    
    cur_edges_1 <- paste(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_edges_2 <- paste(to(colPair(cur_sce)), from(colPair(cur_sce)))
    
    expect_equal(sort(cur_edges_1), sort(cur_edges_2))
    
    expect_equal(sum(isRedundantHit(colPair(cur_sce))), length(colPair(cur_sce))/2)
    
    cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "delaunay", directed = FALSE)
    expect_equal(length(colPair(cur_sce)), 2082)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                              colPairName = "delaunay_interaction_graph")
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4))
    expect_equal(to(colPair(cur_sce))[10:20], c(6, 11, 13, 50,  4,  6, 10, 12,  1,  2,  3))
    
    cur_edges_1 <- paste(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_edges_2 <- paste(to(colPair(cur_sce)), from(colPair(cur_sce)))
    
    expect_equal(sort(cur_edges_1), sort(cur_edges_2))
    
    expect_equal(sum(isRedundantHit(colPair(cur_sce))), length(colPair(cur_sce))/2)
    
    # Max dist
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                               type = "delaunay", max_dist = 20))
    expect_equal(colPairNames(cur_sce), "delaunay_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 1956)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "delaunay_interaction_graph")
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(3L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L))
    expect_equal(to(colPair(cur_sce))[10:20], c(12L, 3L, 7L, 8L, 10L, 2L, 9L, 13L, 14L, 19L, 24L))
    
    cur_edges_1 <- paste(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_edges_2 <- paste(to(colPair(cur_sce)), from(colPair(cur_sce)))
    
    expect_equal(sort(cur_edges_1), sort(cur_edges_2))
    
    expect_equal(sum(isRedundantHit(colPair(cur_sce))), length(colPair(cur_sce))/2)
    
    # KNN
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 5))
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
    
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 5, directed = FALSE))
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
    
    cur_edges_1 <- paste(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_edges_2 <- paste(to(colPair(cur_sce)), from(colPair(cur_sce)))
    
    expect_equal(sort(cur_edges_1), sort(cur_edges_2))
    
    expect_equal(sum(isRedundantHit(colPair(cur_sce))), length(colPair(cur_sce))/2)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as.directed(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as.directed(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as.directed(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    # Other ks
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 10))
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
    
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 10, directed = FALSE))
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
    
    cur_edges_1 <- paste(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_edges_2 <- paste(to(colPair(cur_sce)), from(colPair(cur_sce)))
    
    expect_equal(sort(cur_edges_1), sort(cur_edges_2))
    
    expect_equal(sum(isRedundantHit(colPair(cur_sce))), length(colPair(cur_sce))/2)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as.directed(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as.directed(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:11])})
    cur_ind_real <- graph_from_adj_list(as.list(as.data.frame(cur_ind_real)))
    cur_ind_real <- as.undirected(cur_ind_real)
    cur_ind_real <- as.directed(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    # Other algorithm
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "knn", k = 10))
    expect_silent(cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "knn", k = 10, name = "vptree",
                                 BNPARAM = BiocNeighbors::VptreeParam()))
    expect_silent(cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "knn", k = 10, name = "annoy",
                                 BNPARAM = BiocNeighbors::AnnoyParam()))
    expect_silent(cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "knn", k = 10, name = "hnsw",
                                 BNPARAM = BiocNeighbors::HnswParam()))
    
    expect_equal(colPairNames(cur_sce), 
                 c("knn_interaction_graph", "vptree", "annoy", "hnsw"))
    expect_equal(colPair(cur_sce, "knn_interaction_graph"),
                 colPair(cur_sce, "vptree"))
    
    # Max dist
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                               type = "knn", k = 5, k_max_dist = 10))
    expect_equal(colPairNames(cur_sce), "knn_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 711)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "knn_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    expect_equal(from(colPair(cur_sce))[10:20], c(9, 9, 10, 10, 11, 12, 12, 12, 13, 13, 13))
    expect_equal(to(colPair(cur_sce))[10:20], c(6, 14,  4,  7, 13,  3, 17, 18,  5, 11, 19))
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_dist_real <- apply(as.matrix(cur_dist), 1, function(x){x[sort(order(x)[2:6])]})
    cur_ind_real <- lapply(1:ncol(cur_ind_real), function(x){
        cur_ind_real[cur_dist_real[,x] <= 10,x]
    })
    cur_ind_real <- cur_ind_real[lapply(cur_ind_real, length) > 0]
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){as.numeric(names(which(x == 1)))})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_dist_real <- apply(as.matrix(cur_dist), 1, function(x){x[sort(order(x)[2:6])]})
    cur_ind_real <- lapply(1:ncol(cur_ind_real), function(x){
        cur_ind_real[cur_dist_real[,x] <= 10,x]
    })
    cur_ind_real <- cur_ind_real[lapply(cur_ind_real, length) > 0]
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){as.numeric(names(which(x == 1)))})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
   
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){sort(order(x)[2:6])})
    cur_dist_real <- apply(as.matrix(cur_dist), 1, function(x){x[sort(order(x)[2:6])]})
    cur_ind_real <- lapply(1:ncol(cur_ind_real), function(x){
        cur_ind_real[cur_dist_real[,x] <= 10,x]
    })
    cur_ind_real <- cur_ind_real[lapply(cur_ind_real, length) > 0]
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    cur_ind_test <- apply(cur_ind_test, 1, function(x){as.numeric(names(which(x == 1)))})
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    # Expansion
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "expansion", threshold = 15))
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
    
    cur_edges_1 <- paste(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_edges_2 <- paste(to(colPair(cur_sce)), from(colPair(cur_sce)))
    
    expect_equal(sort(cur_edges_1), sort(cur_edges_2))
    
    expect_equal(sum(isRedundantHit(colPair(cur_sce))), length(colPair(cur_sce))/2)
    
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                 type = "expansion", threshold = 15, directed = FALSE))
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
    
    cur_edges_1 <- paste(from(colPair(cur_sce)), to(colPair(cur_sce)))
    cur_edges_2 <- paste(to(colPair(cur_sce)), from(colPair(cur_sce)))
    
    expect_equal(sort(cur_edges_1), sort(cur_edges_2))
    
    expect_equal(sum(isRedundantHit(colPair(cur_sce))), length(colPair(cur_sce))/2)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 1]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){which(x <= 15)})
    cur_ind_real <- graph_from_adj_list(cur_ind_real)
    cur_ind_real <- simplify(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 2]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){which(x <= 15)})
    cur_ind_real <- graph_from_adj_list(cur_ind_real)
    cur_ind_real <- simplify(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    cur_sce_1 <- cur_sce[,cur_sce$ImageNb == 3]
    
    cur_dist <- dist(as.matrix(colData(cur_sce_1)[,c("Pos_X", "Pos_Y")]))
    cur_ind_real <- apply(as.matrix(cur_dist), 1, function(x){which(x <= 15)})
    cur_ind_real <- graph_from_adj_list(cur_ind_real)
    cur_ind_real <- simplify(cur_ind_real)
    cur_ind_real <- as_edgelist(cur_ind_real)
    cur_ind_real <- table(cur_ind_real[,1], cur_ind_real[,2])
    cur_ind_test <- table(from(colPair(cur_sce_1)), to(colPair(cur_sce_1)))
    
    expect_equal(cur_ind_real, cur_ind_test, check.attributes = FALSE)
    
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                               type = "expansion", threshold = 1))
    expect_equal(colPairNames(cur_sce), "expansion_interaction_graph")
    expect_true(!is.null(colPair(cur_sce)))
    expect_equal(length(colPair(cur_sce)), 0)
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "expansion_interaction_graph")
    expect_silent(print(p))
    p <- plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE,
                     colPairName = "expansion_interaction_graph", 
                     directed = FALSE)
    expect_silent(print(p))
    
    # Parallelisation
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                               type = "expansion", threshold = 15,
                                               BPPARAM = BiocParallel::bpparam()))
    expect_equal(colPairNames(cur_sce), "expansion_interaction_graph")
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                               type = "delaunay",
                                               BPPARAM = BiocParallel::bpparam()))
    expect_equal(colPairNames(cur_sce), "delaunay_interaction_graph")
    expect_silent(cur_sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                               type = "knn", k = 5,
                                               BPPARAM = BiocParallel::bpparam()))
    expect_equal(colPairNames(cur_sce), "knn_interaction_graph")
    
    # SpatialExperiment
    library(SpatialExperiment)
    pancreasSPE <- SpatialExperiment(assays = list(counts = counts(pancreasSCE)),
                                 sample_id = pancreasSCE$ImageName) 
    colData(pancreasSPE) <- colData(pancreasSCE)
    colPairs(pancreasSPE) <- colPairs(pancreasSCE)
    
    expect_error(cur_sce <- buildSpatialGraph(pancreasSPE, img_id = "ImageNb", 
                                               type = "expansion", threshold = 15),
                 regexp = "'coords' not in spatialCoords(object).",
                 fixed = TRUE)
    
    spatialCoords(pancreasSPE) <- as.matrix(colData(pancreasSCE)[,c("Pos_X", "Pos_Y")])
    
    expect_silent(cur_spe <- buildSpatialGraph(pancreasSPE, img_id = "ImageNb", 
                                               type = "expansion", threshold = 15))
    
    expect_equal(colPairNames(cur_spe), "expansion_interaction_graph")
    expect_silent(cur_spe <- buildSpatialGraph(pancreasSPE, img_id = "ImageNb", 
                                               type = "delaunay"))
    expect_equal(colPairNames(cur_spe), "delaunay_interaction_graph")
    expect_silent(cur_spe <- buildSpatialGraph(pancreasSPE, img_id = "ImageNb", 
                                               type = "knn", k = 5))
    expect_equal(colPairNames(cur_spe), "knn_interaction_graph")
    
    # Fail
    expect_error(buildSpatialGraph("test"),
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE),
                 regexp = "argument \"img_id\" is missing, with no default",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = 1),
                 regexp = "'img_id' must be a single string.",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "test"),
                 regexp = "'img_id' not in colData(object).",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "expansion"),
                 regexp = "When constructing a graph via expansion, please specify 'threshold'.",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "expansion"),
                 regexp = "When constructing a graph via expansion, please specify 'threshold'.",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "expansion",
                                   threshold = "test"),
                 regexp = "'threshold' must be a single numeric",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn"),
                 regexp = "When constructing a graph via nearest neighbour detection, please specify 'k'.",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                   k = "test"),
                 regexp = "'k' must be a single numeric",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                   k = 10, k_max_dist = "test"),
                 regexp = "'k_max_dist' must be a single numeric",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "delaunay",
                                   coords = "test"),
                 regexp = "'coords' must be a character vector of length 2.",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "delaunay",
                                   coords = c(2, 1)),
                 regexp = "'coords' must be a character vector of length 2.",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "delaunay",
                                   coords = c("Pos_X", "test")),
                 regexp = "'coords' not in colData(object).",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSPE, img_id = "ImageNb", type = "delaunay",
                                   coords = c("Pos_X", "test")),
                 regexp = "'coords' not in spatialCoords(object).",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "delaunay",
                                   name = 1),
                 regexp = "'name' must be a single string.",
                 fixed = TRUE)
    expect_error(buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "delaunay",
                                   directed = 1),
                 regexp = "'directed' must be a single logical.",
                 fixed = TRUE)
})
