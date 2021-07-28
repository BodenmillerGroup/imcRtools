test_that("aggregateNeighbors function works", {
    library(cytomapper)
    data("pancreasSCE")
    
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "knn",k = 10,
                                     name = "knn_10")
    
    # Works metadata version on knn
    # Directed
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "metadata",
                                                count_by = "CellType"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    expect_true(all(rowSums(as.matrix(cur_sce$aggregatedNeighbors)) == 1))
    
    # check for correct results of neighboring metadata
    cur_dat <- data.frame(colPair(pancreasSCE,"knn_10"))
    cur_dat$celltype <- factor(colData(cur_sce)$CellType[cur_dat$to])
    
    cur_dat <- cur_dat %>% group_by(from) %>% count(celltype, .drop = FALSE) %>% 
        pivot_wider(names_from = "celltype", values_from = "n") %>% ungroup() %>% 
        select(-from) %>% as.matrix()
    
    cur_dat <- cur_dat / rowSums(cur_dat)
    
    expect_equal(cur_sce$aggregatedNeighbors,
                 DataFrame(cur_dat))
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "metadata",
                                                count_by = "CellType",
                                                proportions = FALSE))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    expect_true(all(rowSums(as.matrix(cur_sce$aggregatedNeighbors)) == 10))
    
    # check for correct results of neighboring metadata
    cur_dat <- data.frame(colPair(pancreasSCE,"knn_10"))
    cur_dat$celltype <- factor(colData(cur_sce)$CellType[cur_dat$to])
    
    cur_dat <- cur_dat %>% group_by(from) %>% count(celltype, .drop = FALSE) %>% 
        pivot_wider(names_from = "celltype", values_from = "n") %>% ungroup() %>% 
        select(-from) %>% as.matrix()
    
    expect_equal(cur_sce$aggregatedNeighbors,
                 DataFrame(cur_dat))
    
    cur_sce2 <- cur_sce
    colPair(cur_sce2, "knn_10") <- colPair(cur_sce2, "knn_10")[from(colPair(cur_sce2, "knn_10")) == 215,]
    cur_sce2 <- cur_sce2[,c(215, to(colPair(cur_sce2, "knn_10"))[from(colPair(cur_sce2, "knn_10")) == 215])]
    
    plotSpatial(cur_sce2, img_id = "ImageNb", node_color_by = "CellType", draw_edges = TRUE,
                colPairName = "knn_10")
    cur_sce2$aggregatedNeighbors
    
    cur_sce2 <- cur_sce
    colPair(cur_sce2, "knn_10") <- colPair(cur_sce2, "knn_10")[from(colPair(cur_sce2, "knn_10")) == 10,]
    cur_sce2 <- cur_sce2[,c(10, to(colPair(cur_sce2, "knn_10"))[from(colPair(cur_sce2, "knn_10")) == 10])]
    
    plotSpatial(cur_sce2, img_id = "ImageNb", node_color_by = "CellType", draw_edges = TRUE,
                colPairName = "knn_10")
    cur_sce2$aggregatedNeighbors
    
    cur_sce2 <- cur_sce
    colPair(cur_sce2, "knn_10") <- colPair(cur_sce2, "knn_10")[from(colPair(cur_sce2, "knn_10")) == 322,]
    cur_sce2 <- cur_sce2[,c(322, to(colPair(cur_sce2, "knn_10"))[from(colPair(cur_sce2, "knn_10")) == 322])]
    
    plotSpatial(cur_sce2, img_id = "ImageNb", node_color_by = "CellType", draw_edges = TRUE,
                colPairName = "knn_10")
    cur_sce2$aggregatedNeighbors
    
    # Undirected
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "knn",k = 10,
                                     directed = FALSE,
                                     name = "knn_10")
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "metadata",
                                                count_by = "CellType"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    expect_true(all(rowSums(as.matrix(cur_sce$aggregatedNeighbors)) == 1))
    
    # check for correct results of neighboring metadata
    cur_dat <- data.frame(colPair(pancreasSCE,"knn_10"))
    cur_dat$celltype <- factor(colData(cur_sce)$CellType[cur_dat$to])
    
    cur_dat <- cur_dat %>% group_by(from) %>% count(celltype, .drop = FALSE) %>% 
        pivot_wider(names_from = "celltype", values_from = "n") %>% ungroup() %>% 
        select(-from) %>% as.matrix()
    
    cur_dat <- cur_dat / rowSums(cur_dat)
    
    expect_equal(cur_sce$aggregatedNeighbors,
                 DataFrame(cur_dat))
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "metadata",
                                                count_by = "CellType",
                                                proportions = FALSE))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    expect_true(all(rowSums(as.matrix(cur_sce$aggregatedNeighbors)) == countLnodeHits(colPair(cur_sce, "knn_10"))))
    
    # check for correct results of neighboring metadata
    cur_dat <- data.frame(colPair(pancreasSCE,"knn_10"))
    cur_dat$celltype <- factor(colData(cur_sce)$CellType[cur_dat$to])
    
    cur_dat <- cur_dat %>% group_by(from) %>% count(celltype, .drop = FALSE) %>% 
        pivot_wider(names_from = "celltype", values_from = "n") %>% ungroup() %>% 
        select(-from) %>% as.matrix()
    
    expect_equal(cur_sce$aggregatedNeighbors,
                 DataFrame(cur_dat))
    
    cur_sce2 <- cur_sce
    colPair(cur_sce2, "knn_10") <- colPair(cur_sce2, "knn_10")[from(colPair(cur_sce2, "knn_10")) == 215,]
    cur_sce2 <- cur_sce2[,c(215, to(colPair(cur_sce2, "knn_10"))[from(colPair(cur_sce2, "knn_10")) == 215])]
    
    plotSpatial(cur_sce2, img_id = "ImageNb", node_color_by = "CellType", draw_edges = TRUE,
                colPairName = "knn_10")
    cur_sce2$aggregatedNeighbors
    
    cur_sce2 <- cur_sce
    colPair(cur_sce2, "knn_10") <- colPair(cur_sce2, "knn_10")[from(colPair(cur_sce2, "knn_10")) == 10,]
    cur_sce2 <- cur_sce2[,c(10, to(colPair(cur_sce2, "knn_10"))[from(colPair(cur_sce2, "knn_10")) == 10])]
    
    plotSpatial(cur_sce2, img_id = "ImageNb", node_color_by = "CellType", draw_edges = TRUE,
                colPairName = "knn_10")
    cur_sce2$aggregatedNeighbors
    
    cur_sce2 <- cur_sce
    colPair(cur_sce2, "knn_10") <- colPair(cur_sce2, "knn_10")[from(colPair(cur_sce2, "knn_10")) == 322,]
    cur_sce2 <- cur_sce2[,c(322, to(colPair(cur_sce2, "knn_10"))[from(colPair(cur_sce2, "knn_10")) == 322])]
    
    plotSpatial(cur_sce2, img_id = "ImageNb", node_color_by = "CellType", draw_edges = TRUE,
                colPairName = "knn_10")
    cur_sce2$aggregatedNeighbors
    
    # Works expression version on knn
    # Directed
    # Mean
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "knn",k = 10,
                                     name = "knn_10")
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "exprs"))
    
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "mean_aggregatedExpression"))
    expect_s4_class(cur_sce$mean_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$mean_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"knn_10"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"knn_10"))

    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "mean", 
                                      use.assay.type = "exprs")
    
    expect_equal(as.matrix(cur_sce$mean_aggregatedExpression),
                 t(as.matrix(assay(cur_sce_2, "exprs"))), check.attributes = FALSE)
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "counts"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "mean_aggregatedExpression"))
    expect_s4_class(cur_sce$mean_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$mean_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"knn_10"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"knn_10"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "mean", 
                                      use.assay.type = "counts")
    
    expect_equal(as.matrix(cur_sce$mean_aggregatedExpression),
                 t(as.matrix(assay(cur_sce_2, "counts"))), check.attributes = FALSE)
    
    # Median
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "knn",k = 10,
                                     name = "knn_10")
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "exprs",
                                                statistic = "median"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "median_aggregatedExpression"))
    expect_s4_class(cur_sce$median_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$median_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"knn_10"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"knn_10"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "median", 
                                      use.assay.type = "exprs")
    
    expect_equal(as.matrix(cur_sce$median_aggregatedExpression),
                 t(as.matrix(assay(cur_sce_2, "exprs"))), check.attributes = FALSE)
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "counts",
                                                statistic = "median"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "median_aggregatedExpression"))
    expect_s4_class(cur_sce$median_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$median_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"knn_10"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"knn_10"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "median", 
                                      use.assay.type = "counts")
    
    expect_equal(as.matrix(cur_sce$median_aggregatedExpression),
                 t(as.matrix(assay(cur_sce_2, "counts"))), check.attributes = FALSE)
    
    # Undirected
    ## Mean
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "knn",k = 10,
                                     name = "knn_10", directed = FALSE)
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "exprs"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "mean_aggregatedExpression"))
    expect_s4_class(cur_sce$mean_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$mean_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"knn_10"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"knn_10"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "mean", 
                                      use.assay.type = "exprs")
    
    expect_equal(as.matrix(cur_sce$mean_aggregatedExpression),
                 t(as.matrix(assay(cur_sce_2, "exprs"))), check.attributes = FALSE)
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "counts"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "mean_aggregatedExpression"))
    expect_s4_class(cur_sce$mean_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$mean_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"knn_10"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"knn_10"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "mean", 
                                      use.assay.type = "counts")
    
    expect_equal(as.matrix(cur_sce$mean_aggregatedExpression),
                 t(as.matrix(assay(cur_sce_2, "counts"))), check.attributes = FALSE)
    
    # Median
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "knn",k = 10,
                                     name = "knn_10",
                                     directed = FALSE)
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "exprs",
                                                statistic = "median"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "median_aggregatedExpression"))
    expect_s4_class(cur_sce$median_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$median_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"knn_10"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"knn_10"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "median", 
                                      use.assay.type = "exprs")
    
    expect_equal(as.matrix(cur_sce$median_aggregatedExpression),
                 t(as.matrix(assay(cur_sce_2, "exprs"))), check.attributes = FALSE)
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "counts",
                                                statistic = "median"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "median_aggregatedExpression"))
    expect_s4_class(cur_sce$median_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$median_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"knn_10"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"knn_10"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "median", 
                                      use.assay.type = "counts")
    
    expect_equal(as.matrix(cur_sce$median_aggregatedExpression),
                 t(as.matrix(assay(cur_sce_2, "counts"))), check.attributes = FALSE)
    
    # Expansion
    ## metadata
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "expansion",
                                     threshold = 5,
                                     name = "exp_20")
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    
    # check for correct results of neighboring metadata
    cur_dat <- data.frame(colPair(pancreasSCE,"exp_20"))
    cur_dat$celltype <- factor(colData(cur_sce)$CellType)[cur_dat$to]
    
    cur_dat <- cur_dat %>% group_by(from) %>% count(celltype, .drop = FALSE) %>% 
        pivot_wider(names_from = "celltype", values_from = "n") %>% ungroup() %>% as.matrix()
    
    cur_dat[,-1] <- cur_dat[,-1] / rowSums(cur_dat[,-1])
    
    expect_equal(cur_sce$aggregatedNeighbors[cur_dat[,"from"],],
                 DataFrame(cur_dat[,-1]))
    expect_true(all(as.matrix(cur_sce$aggregatedNeighbors[-cur_dat[,"from"],]) == 0))
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType",
                                                proportions = FALSE))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    expect_true(all(rowSums(as.matrix(cur_sce$aggregatedNeighbors)) == countLnodeHits(colPair(cur_sce, "exp_20"))))
    
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "expansion",
                                     threshold = 10,
                                     name = "exp_20")
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    
    # check for correct results of neighboring metadata
    cur_dat <- data.frame(colPair(pancreasSCE,"exp_20"))
    cur_dat$celltype <- factor(colData(cur_sce)$CellType)[cur_dat$to]
    
    cur_dat <- cur_dat %>% group_by(from) %>% count(celltype, .drop = FALSE) %>% 
        pivot_wider(names_from = "celltype", values_from = "n") %>% ungroup() %>% as.matrix()
    
    cur_dat[,-1] <- cur_dat[,-1] / rowSums(cur_dat[,-1])
    
    expect_equal(cur_sce$aggregatedNeighbors[cur_dat[,"from"],],
                 DataFrame(cur_dat[,-1]))
    expect_true(all(as.matrix(cur_sce$aggregatedNeighbors[-cur_dat[,"from"],]) == 0))
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType",
                                                proportions = FALSE))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    expect_true(all(rowSums(as.matrix(cur_sce$aggregatedNeighbors)) == countLnodeHits(colPair(cur_sce, "exp_20"))))
    
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "expansion",
                                     threshold = 15,
                                     name = "exp_20")
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    
    # check for correct results of neighboring metadata
    cur_dat <- data.frame(colPair(pancreasSCE,"exp_20"))
    cur_dat$celltype <- factor(colData(cur_sce)$CellType)[cur_dat$to]
    
    cur_dat <- cur_dat %>% group_by(from) %>% count(celltype, .drop = FALSE) %>% 
        pivot_wider(names_from = "celltype", values_from = "n") %>% ungroup() %>% as.matrix()
    
    cur_dat[,-1] <- cur_dat[,-1] / rowSums(cur_dat[,-1])
    
    expect_equal(cur_sce$aggregatedNeighbors[cur_dat[,"from"],],
                 DataFrame(cur_dat[,-1]))
    expect_true(all(as.matrix(cur_sce$aggregatedNeighbors[-cur_dat[,"from"],]) == 0))
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType",
                                                proportions = FALSE))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    expect_true(all(rowSums(as.matrix(cur_sce$aggregatedNeighbors)) == countLnodeHits(colPair(cur_sce, "exp_20"))))
    
    ## expression
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "expansion",
                                     threshold = 5,
                                     name = "exp_20")
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "expression",
                                                assay_type = "exprs"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "mean_aggregatedExpression"))
    expect_s4_class(cur_sce$mean_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$mean_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"exp_20"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"exp_20"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "mean", 
                                      use.assay.type = "exprs")
    
    expect_equal(as.matrix(cur_sce$mean_aggregatedExpression[unique(from(colPair(pancreasSCE,"exp_20"))),]),
                 t(as.matrix(assay(cur_sce_2, "exprs"))), check.attributes = FALSE)
    expect_true(all(is.na(as.matrix(cur_sce$mean_aggregatedExpression[-unique(from(colPair(pancreasSCE,"exp_20"))),]))))
    
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "expression",
                                                assay_type = "counts"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "mean_aggregatedExpression"))
    expect_s4_class(cur_sce$mean_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$mean_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"exp_20"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"exp_20"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "mean", 
                                      use.assay.type = "counts")
    
    expect_equal(as.matrix(cur_sce$mean_aggregatedExpression[from(colPair(pancreasSCE,"exp_20")),]),
                 t(as.matrix(assay(cur_sce_2, "counts"))), check.attributes = FALSE)
    expect_true(all(is.na(as.matrix(cur_sce$mean_aggregatedExpression[-from(colPair(pancreasSCE,"exp_20")),]))))
    
    # Median
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "expression",
                                                assay_type = "exprs",
                                                statistic = "median"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "median_aggregatedExpression"))
    expect_s4_class(cur_sce$median_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$median_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"exp_20"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"exp_20"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "median", 
                                      use.assay.type = "exprs")
    
    expect_equal(as.matrix(cur_sce$median_aggregatedExpression[unique(from(colPair(pancreasSCE,"exp_20"))),]),
                 t(as.matrix(assay(cur_sce_2, "exprs"))), check.attributes = FALSE)
    expect_true(all(is.na(as.matrix(cur_sce$median_aggregatedExpression[-unique(from(colPair(pancreasSCE,"exp_20"))),]))))
    
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "expansion",
                                     threshold = 15,
                                     name = "exp_20")
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "expression",
                                                assay_type = "exprs"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "mean_aggregatedExpression"))
    expect_s4_class(cur_sce$mean_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$mean_aggregatedExpression)), dim(cur_sce))
    
    # check for correct results of neighboring metadata
    cur_sce_2 <- cur_sce[,to(colPair(pancreasSCE,"exp_20"))]
    colData(cur_sce_2)$from <- from(colPair(pancreasSCE,"exp_20"))
    
    cur_sce_2 <- aggregateAcrossCells(cur_sce_2, ids = cur_sce_2$from, statistics = "mean", 
                                      use.assay.type = "exprs")
    
    expect_equal(as.matrix(cur_sce$mean_aggregatedExpression[unique(from(colPair(pancreasSCE,"exp_20"))),]),
                 t(as.matrix(assay(cur_sce_2, "exprs"))), check.attributes = FALSE)
    expect_true(all(is.na(as.matrix(cur_sce$mean_aggregatedExpression[-unique(from(colPair(pancreasSCE,"exp_20"))),]))))
    
    
    # Delaunay
    ## metadata
    
    ## expression
    
    
    
    # Works metadata version on expansion
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType"))
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType",
                                                name = "exp_20_nb"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    
    # Works expression version on expansion
    cur_markers <- rownames(pancreasSCE)
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "expression",
                                                assay_type = "exprs",
                                                subset_row = cur_markers,
                                                name = "exp_20_expMean"))
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "expression",
                                                assay_type = "exprs",
                                                subset_row = cur_markers,
                                                name = "exp_20_expMedian",
                                                statistic = "median"))
    
   
    # Error
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "expansion",
                                     threshold = 1,
                                     name = "exp_20")
    expect_error(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "exp_20",
                                                aggregate_by = "metadata",
                                                count_by = "CellType"),
                 regexp = "No interactions found.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors("test"),
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE),
                 regexp = "argument \"colPairName\" is missing, with no default",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName =  "test"),
                 regexp = "'colPairName' not in 'colPairNames(object)'.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "metadata"),
                 regexp = "Provide a 'colData(object)' entry to aggregate by.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "metadata", count_by = "test"),
                 regexp = "'count_by' is not a valid enty of 'colData(object)'.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "expression"),
                 regexp = "'assay_type' not provided",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "expression", assay_type = "test"),
                 regexp = "'assay_type' not an assay in the 'object'.",
                 fixed = TRUE)
  
})
