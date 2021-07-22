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
    
    # Works celltypes version on knn
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "knn_10",summarize_by = "celltypes",
                                                    group = "CellType",name = "knn_10_nb"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    
    # Works expression version on knn
    cur_markers <- rownames(pancreasSCE)
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "knn_10",summarize_by = "expression",
                                                    assay_type = "exprs",subset_row = cur_markers,name = "knn_10_expMean"))
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "knn_10",summarize_by = "expression",
                                                    assay_type = "exprs",subset_row = cur_markers,name = "knn_10_expMedian",summaryStats = "median"))
    
    # Works celltypes version on expansion
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "exp_20",summarize_by = "celltypes",
                                                    group = "CellType"))
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "exp_20",summarize_by = "celltypes",
                                                    group = "CellType",name = "exp_20_nb"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    
    # Works expression version on expansion
    cur_markers <- rownames(pancreasSCE)
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "exp_20",summarize_by = "expression",
                                                    assay_type = "exprs",subset_row = cur_markers,name = "exp_20_expMean"))
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "exp_20",summarize_by = "expression",
                                                    assay_type = "exprs",subset_row = cur_markers,name = "exp_20_expMedian",summaryStats = "median"))
    
    # check for correct results of neighboring celltypes
    cur_dat <- data.frame(colPair(pancreasSCE,"knn_10"))
    # our reference cell will be cell number 215
    cur_dat <- cur_dat[which(cur_dat$from == 215),]
    # what type of cells are the neighbors of cell 71?
    neighbors_cell_215 <- unclass(as.matrix(t(table(pancreasSCE[,cur_dat$to]$CellType))))

    
    cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "knn_10",summarize_by = "celltypes",
                       group = "CellType")
    
    expect_equal(as.data.frame(cur_sce$summarizedNeighbors[215,]),as.data.frame(neighbors_cell_215))
    
    
    # check for correct results of neighboring expression
    cur_dat <- data.frame(colPair(pancreasSCE,"knn_10"))
    # our reference cell will be cell number 215
    cur_dat <- cur_dat[which(cur_dat$from == 215),]
    # what type of cells are the neighbors of cell 215?
    neighbors_cell_215 <- t(assay(pancreasSCE,"exprs"))[cur_dat$to,]
    # calculate mean per marker for neighboring cells
    mean_exp <- apply(neighbors_cell_215,2,mean)

    cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "knn_10",summarize_by = "expression",
                                  assay_type = "exprs",subset_row = cur_markers,name = "knn_10_expMean")
    
    expect_equal(as.numeric(as.data.frame(cur_sce$knn_10_expMean[215,])),as.numeric(mean_exp))
    
    # calculate median per marker for neighboring cells
    median_exp <- apply(neighbors_cell_215,2,median)
    
    cur_sce <- aggregateNeighbors(object = pancreasSCE,colPairName = "knn_10",summarize_by = "expression",
                                  assay_type = "exprs",subset_row = cur_markers,name = "knn_10_expMean",summaryStats = "median")
    
    expect_equal(as.numeric(as.data.frame(cur_sce$knn_10_expMean[215,])),as.numeric(median_exp))
    
    # Error
    expect_error(aggregateNeighbors("test"),
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE),
                 regexp = "argument \"colPairName\" is missing, with no default",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName =  "test"),
                 regexp = "'colPairName' not in colPair(object).",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",summarize_by = "celltypes"),
                 regexp = "provide a colData entry to summarize by",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",summarize_by = "celltypes", group = "test"),
                 regexp = "'group' is not a valid enty of colData(object).",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",summarize_by = "expression"),
                 regexp = "'assay_type' not provided",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",summarize_by = "expression", assay_type = "test"),
                 regexp = "'assay_type' not an assay in the 'object'.",
                 fixed = TRUE)
  
})

#blabbla