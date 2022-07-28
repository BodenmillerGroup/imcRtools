test_that("patchDetection function works", {
    library(cytomapper)
    data(pancreasSCE)

    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                    type = "expansion", threshold = 20)

    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                   patch_cells = pancreasSCE$CellType == "celltype_B",
                                   colPairName = "expansion_interaction_graph"))
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph"))
    
    # SpatialExperiment
    cur_spe <- SpatialExperiment:::.sce_to_spe(cur_sce, sample_id = as.character(pancreasSCE$ImageNb))
    spatialCoords(cur_spe) <- as.matrix(colData(cur_sce)[,c("Pos_X", "Pos_Y")])
    colData(cur_spe)[c("Pos_X", "Pos_Y")] <- NULL
    
    
    
    
    # Check for sparse graphs
    cur_sce <- pancreasSCE
    cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                     type = "expansion", threshold = 10)
    
    expect_silent(cur_sce <- patchDetection(cur_sce, 
                                            patch_cells = cur_sce$CellType == "celltype_C",
                                            colPairName = "expansion_interaction_graph"))
    
    # Error
    expect_error(patchDetection("test"),
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = "test"),
                 regexp = "'patch_cells' must all be logical.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = TRUE),
                 regexp = "Length of 'patch_cells' must match the number of cells in 'object'.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = 1),
                 regexp = "'colPairName' must be a single string.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "test"),
                 regexp = "'colPairName' not in 'colPairNames(object)'.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                min_patch_size = c(1,2)),
                 regexp = "'min_patch_size' must be a single numeric.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                coords = 1),
                 regexp = "'coords' must be a character vector of length 2.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                coords = c("test", "Pos_y")),
                 regexp = "'coords' not in colData(object).",
                 fixed = TRUE)
    expect_error(patchDetection(cur_spe, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                coords = c("test", "Pos_y")),
                 regexp = "'coords' not in spatialCoords(object).",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                name = c("test", "Pos_y")),
                 regexp = "'name' must be a single string.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                expand_by = "test"),
                 regexp = "'expand_by' must be a single numeric.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                convex = "test"),
                 regexp = "'convex' must be a single logical.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                expand_by = 10),
                 regexp = "'img_id' must be specified when patch expansion is performed.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                expand_by = 10, img_id = 1),
                 regexp = "'img_id' must be a single string.",
                 fixed = TRUE)
    expect_error(patchDetection(pancreasSCE, patch_cells = pancreasSCE$Pattern,
                                colPairName = "expansion_interaction_graph",
                                expand_by = 10, img_id = "test"),
                 regexp = "'img_id' not in colData(object).",
                 fixed = TRUE)
                     
})
