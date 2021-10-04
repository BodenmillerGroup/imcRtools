test_that("aggregateNeighbors function works", {
    library(cytomapper)
    library(dplyr)
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
    
    # Delaunay
    ## metadata
    pancreasSCE <- buildSpatialGraph(object = pancreasSCE,
                                     img_id = "ImageNb",
                                     type = "delaunay")
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "delaunay_interaction_graph",
                                                aggregate_by = "metadata",
                                                count_by = "CellType"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "aggregatedNeighbors"))
    expect_s4_class(cur_sce$aggregatedNeighbors, "DataFrame")
    
    ## expression
    expect_silent(cur_sce <- aggregateNeighbors(object = pancreasSCE,
                                                colPairName = "delaunay_interaction_graph",
                                                aggregate_by = "expression",
                                                assay_type = "exprs"))
    
    expect_s4_class(cur_sce , class = "SingleCellExperiment")
    expect_equal(names(colData(cur_sce)), 
                 c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb",
                   "CellNb", "MaskName", "Pattern", "mean_aggregatedExpression"))
    expect_s4_class(cur_sce$mean_aggregatedExpression, "DataFrame")
    expect_equal(dim(t(cur_sce$mean_aggregatedExpression)), dim(cur_sce))
    
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
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName =  c("test", 1)),
                 regexp = "'colPairName' must be a single string.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName =  "test"),
                 regexp = "'colPairName' not in 'colPairNames(object)'.",
                 fixed = TRUE)
    
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "metadata"),
                 regexp = "Provide a 'colData(object)' entry to aggregate by.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "metadata", count_by = 1),
                 regexp = "'count_by' must be a single string.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "metadata", count_by = c("CellType", "test")),
                 regexp = "'count_by' must be a single string.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "metadata", count_by = "test"),
                 regexp = "'count_by' is not a valid enty of 'colData(object)'.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "metadata", count_by = "CellType",
                                    proportions = "test"),
                 regexp = "'proportions' must be a single logical",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "expression"),
                 regexp = "'assay_type' not provided",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "expression", assay_type = 1),
                 regexp = "'assay_type' must be a single string.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "expression", assay_type = "test"),
                 regexp = "'assay_type' not an assay in the 'object'.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "expression", assay_type = "counts",
                                    subset_row = c(1, "test")),
                 regexp = "'subset_row' not in rownames(object).",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "expression", assay_type = "counts",
                                    subset_row = TRUE),
                 regexp = "'subset_row' logical entries must be as long as 'nrow(object)'.",
                 fixed = TRUE)
    expect_error(aggregateNeighbors(object = pancreasSCE, colPairName = "knn_10",
                                    aggregate_by = "expression", assay_type = "counts",
                                    name = TRUE),
                 regexp = "'name' must be a single string.",
                 fixed = TRUE)
})
