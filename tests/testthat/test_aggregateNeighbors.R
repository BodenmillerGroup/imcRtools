test_that("aggregateNeighbors function works", {
    library(cytomapper)
    data("pancreasSCE")
    
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "knn",k = 10,
                                     name = "knn_10")
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "expansion",
                                     threshold = 20,
                                     name = "exp_20")
    
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
    
    
    # Works expression version on knn
    cur_markers <- rownames(pancreasSCE)
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "exprs",
                                                subset_row = cur_markers,
                                                name = "knn_10_expMean"))
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "knn_10",
                                                aggregate_by = "expression",
                                                assay_type = "exprs",
                                                subset_row = cur_markers,
                                                name = "knn_10_expMedian",
                                                statistic = "median"))
    
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
    
    # check for correct results of neighboring metadata
    cur_dat <- data.frame(colPair(pancreasSCE,"knn_10"))
    # our reference cell will be cell number 215
    cur_dat <- cur_dat[which(cur_dat$from == 215),]
    # what type of cells are the neighbors of cell 71?
    neighbors_cell_215 <- unclass(as.matrix(t(table(pancreasSCE[,cur_dat$to]$CellType))))

    
    cur_sce <- aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                  aggregate_by = "metadata",
                                  count_by = "CellType")
    
    expect_equal(as.data.frame(cur_sce$aggregatedNeighbors[215,]),
                 as.data.frame(neighbors_cell_215/sum(neighbors_cell_215)))
    
    
    # check for correct results of neighboring expression
    cur_dat <- data.frame(colPair(pancreasSCE,"knn_10"))
    # our reference cell will be cell number 215
    cur_dat <- cur_dat[which(cur_dat$from == 215),]
    # what type of cells are the neighbors of cell 215?
    neighbors_cell_215 <- t(assay(pancreasSCE,"exprs"))[cur_dat$to,]
    # calculate mean per marker for neighboring cells
    mean_exp <- apply(neighbors_cell_215,2,mean)

    cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "knn_10",
                                  aggregate_by = "expression",
                                  assay_type = "exprs",
                                  subset_row = cur_markers,
                                  name = "knn_10_expMean")
    
    expect_equal(as.numeric(as.data.frame(cur_sce$knn_10_expMean[215,])),
                 as.numeric(mean_exp))
    
    # calculate median per marker for neighboring cells
    median_exp <- apply(neighbors_cell_215,2,median)
    
    cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                  colPairName = "knn_10",
                                  aggregate_by = "expression",
                                  assay_type = "exprs",
                                  subset_row = cur_markers,
                                  name = "knn_10_expMean",
                                  statistic = "median")
    
    expect_equal(as.numeric(as.data.frame(cur_sce$knn_10_expMean[215,])),
                 as.numeric(median_exp))
    
    # Error
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
