test_that("patchDetection function works", {
    library(cytomapper)
    data(pancreasSCE)

    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                    type = "expansion", threshold = 20)

    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                   patch_cells = pancreasSCE$CellType == "celltype_B",
                                   colPairName = "expansion_interaction_graph"))
    
    expect_true(is(cur_sce, "SingleCellExperiment"))
    expect_true("patch_id" %in% names(colData(cur_sce)))
    expect_equal(unique(cur_sce$patch_id), c(NA, "1", "2", "3", "4", "5", "6", "7", "8"))
    
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "CellType")  
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            min_patch_size = 5))
    
    expect_equal(unique(cur_sce$patch_id), c(NA, "6", "7", "8"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType %in% c("celltype_B", "celltype_A"),
                                            colPairName = "expansion_interaction_graph"))
    expect_equal(unique(cur_sce$patch_id), c(NA, "1", "2", "3"))
    
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType %in% c("celltype_B", "celltype_A", "celltype_C"),
                                            colPairName = "expansion_interaction_graph"))
    expect_equal(unique(cur_sce$patch_id), c("1", "2", "3"))
    
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 10, img_id = "ImageNb",
                                            name = "patch_id_2"))
    expect_equal(unique(cur_sce$patch_id_2), c(NA, "1", "2", "3", "4", "5", "6", "7", "8"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id_2")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 50, img_id = "ImageNb",
                                            name = "patch_id"))
    expect_equal(unique(cur_sce$patch_id), c(NA,  "1", "3", "6", "7", "8"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 1000, img_id = "ImageNb",
                                            name = "patch_id"))
    expect_equal(unique(cur_sce$patch_id), c("3", "6", "7", "8"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 20, img_id = "ImageNb",
                                            name = "patch_id",
                                            min_patch_size = 5))
    expect_equal(unique(cur_sce$patch_id), c(NA, "6", "7", "8"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    cur_sce_2 <- pancreasSCE
    cur_sce_2$X <- cur_sce_2$Pos_X
    cur_sce_2$Y <- cur_sce_2$Pos_Y    
    cur_sce_2$Pos_X <- NULL
    cur_sce_2$Pos_y <- NULL
    
    expect_silent(cur_sce <- patchDetection(cur_sce_2, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 30,
                                            img_id = "ImageNb",
                                            coords = c("X", "Y")))
    
    expect_equal(unique(cur_sce$patch_id), c(NA, "1", "3", "4", "6", "7", "8"))
    
    # Concave and convex
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 1, img_id = "ImageNb",
                                            name = "patch_id",
                                            convex = TRUE))
    expect_equal(unique(cur_sce$patch_id), c(NA, "1", "2", "3", "4", "5", "6", "7", "8"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 20, img_id = "ImageNb",
                                            name = "patch_id",
                                            convex = TRUE))
    expect_equal(unique(cur_sce$patch_id), c(NA, "1", "2", "3", "4", "5", "6", "7", "8"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 1000, img_id = "ImageNb",
                                            name = "patch_id",
                                            convex = TRUE))
    expect_equal(unique(cur_sce$patch_id), c("3", "6", "7", "8"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    # Other graph constructors
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                     type = "expansion", threshold = 1)
    
    expect_error(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph"),
                 regexp = "No interactions found.",
                 fixed = TRUE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                     type = "expansion", threshold = 5)
    expect_error(cur_sce <- patchDetection(pancreasSCE, 
                                           patch_cells = pancreasSCE$CellType == "celltype_B",
                                           colPairName = "expansion_interaction_graph"),
                 regexp = "No connected components found.",
                 fixed = TRUE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                     type = "expansion", threshold = 10)
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_equal(unique(cur_sce$patch_id), c(NA, as.character(seq_len(39))))
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                     type = "knn", k = 3)
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "knn_interaction_graph"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_equal(unique(cur_sce$patch_id), c(NA, as.character(seq_len(13))))
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                     type = "delaunay")
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "delaunay_interaction_graph"))
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_equal(unique(cur_sce$patch_id), c(NA, as.character(seq_len(9))))
    
    # Spatial Experiment
    cur_spe <- SpatialExperiment:::.sce_to_spe(pancreasSCE, sample_id = as.character(pancreasSCE$ImageNb))
    spatialCoords(cur_spe) <- as.matrix(colData(pancreasSCE)[,c("Pos_X", "Pos_Y")])
    
    cur_spe <- buildSpatialGraph(cur_spe, img_id = "ImageNb", 
                                     type = "expansion", threshold = 20)
    
    expect_silent(cur_spe_2 <- patchDetection(cur_spe, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 10, img_id = "ImageNb",
                                            name = "patch_id_2"))
    expect_equal(unique(cur_spe_2$patch_id_2), c(NA, "1", "2", "3", "4", "5", "6", "7", "8"))
    plotSpatial(cur_spe_2, img_id = "ImageNb", node_color_by = "patch_id_2")
    
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
