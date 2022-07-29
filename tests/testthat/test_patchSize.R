test_that("patchSize function works", {
    library(cytomapper)
    data(pancreasSCE)

    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
                                    type = "expansion", threshold = 20)

    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                   patch_cells = pancreasSCE$CellType == "celltype_B",
                                   colPairName = "expansion_interaction_graph"))
    
    expect_silent(cur_out_1 <- patchSize(cur_sce))
    expect_s4_class(cur_out_1, "DataFrame")
    expect_equal(cur_out_1$size, c(NA, NA, 32.6362638000003, NA, NA, 2028.226120055, 169.684533615001, 
                                   3778.12476426999))
    expect_equal(cur_out_1$polygon[3][[1]], structure(list(structure(c(195.8023, 175.9615, 181.2143, 195.8023, 
                                                                       286.2442, 266.8654, 275.2857, 286.2442), dim = c(4L, 2L), dimnames = list(
                                                                           NULL, c("V1", "V2")))), class = c("XY", "POLYGON", "sfg")))
    
    # SpatialExperiment
    cur_spe <- SpatialExperiment:::.sce_to_spe(cur_sce, sample_id = as.character(cur_sce$ImageNb))
    spatialCoords(cur_spe) <- as.matrix(colData(cur_sce)[,c("Pos_X", "Pos_Y")])
    colData(cur_spe)[c("Pos_X", "Pos_Y")] <- NULL
    
    expect_silent(cur_out_2 <- patchSize(cur_spe))
    expect_s4_class(cur_out_2, "DataFrame")
    
    expect_equal(cur_out_1, cur_out_2)
    
    # Check for sparse graphs
    cur_sce_2 <- pancreasSCE
    cur_sce_2 <- buildSpatialGraph(cur_sce_2, img_id = "ImageNb", 
                                     type = "expansion", threshold = 10)
    
    expect_silent(cur_sce_2 <- patchDetection(cur_sce_2, 
                                            patch_cells = cur_sce$CellType == "celltype_C",
                                            colPairName = "expansion_interaction_graph"))
    
    expect_silent(cur_out <- patchSize(cur_sce_2))
    expect_equal(cur_out$size, c(NA, NA, NA, NA, NA, 14.1314614, NA, NA, 512.378910555, NA, 
                                 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1301.014577695, 
                                 NA, NA, NA, 37.3505263149999, 255.10388471, NA, NA, NA, 21.2789770099999, 
                                 NA, 55.1865428500002, 71.2887023200001, NA, NA, 4.99787407499997, 
                                 NA, NA, 180.80620656, NA, NA, NA, 52.6575400250001, NA, NA, NA, 
                                 NA, NA, 644.063649344998, NA, NA, 186.640773035, NA, 31.75246717, 
                                 NA, NA, 119.468629760002, 208.93155685, NA, NA, NA, NA, NA, NA, 
                                 NA, NA, NA))
    
    # Convex
    expect_silent(cur_out_3 <- patchSize(cur_sce, convex = TRUE))
    expect_false(identical(cur_out_1$size, cur_out_3$size))
    expect_s4_class(cur_out_3, "DataFrame")
    expect_equal(cur_out_3$size, c(NA, NA, 32.6370112445557, NA, NA, 4786.46867623013, 195.870087122285, 
                                   7075.35257151764))
    expect_equal(cur_out_3$polygon[3][[1]], structure(list(structure(c(175.961538461538, 181.214285714286, 
                                                                       195.802325581395, 175.961538461538, 266.865384615385, 275.285714285714, 
                                                                       286.244186046512, 266.865384615385), dim = c(4L, 2L), dimnames = list(
                                                                           NULL, c("Pos_X", "Pos_Y")))), class = c("XY", "POLYGON", 
                                                                                                                   "sfg")))
    
    # Expansion setting
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 10,
                                            img_id = "ImageNb"))
    
    expect_silent(cur_out_4 <- patchSize(cur_sce))
    expect_false(identical(cur_out_1$size, cur_out_3$size))
    expect_s4_class(cur_out_4, "DataFrame")
    expect_equal(cur_out_4$size, c(NA, NA, 95.8792097900002, NA, NA, 2818.688752595, 621.836906870001, 
                                   4742.921943955))
    expect_equal(cur_out_4$polygon[3][[1]], structure(list(structure(c(186.1754, 175.9615, 181.2143, 189.2317, 
                                                                       195.8023, 186.1754, 268.8947, 266.8654, 275.2857, 279.1707, 286.2442, 
                                                                       268.8947), dim = c(6L, 2L), dimnames = list(NULL, c("V1", "V2"
                                                                       )))), class = c("XY", "POLYGON", "sfg")))
    
    # Error
    expect_error(patchSize("test"),
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(patchSize(pancreasSCE, patch_name = "test"),
                 regexp = "'patch_name' nor in 'colData(object)'.",
                 fixed = TRUE)
    expect_error(patchSize(pancreasSCE, patch_name = c("test", "test2")),
                 regexp = "'patch_name' must be a single string.",
                 fixed = TRUE)
    expect_error(patchSize(cur_sce, coords = "test"),
                 regexp = "'coords' must be a character vector of length 2.",
                 fixed = TRUE)
    expect_error(patchSize(cur_sce, coords = c("test1", "test2")),
                 regexp = "'coords' not in colData(object).",
                 fixed = TRUE)
    expect_error(patchSize(cur_spe, coords = c("test1", "test2")),
                 regexp = "'coords' not in spatialCoords(object).",
                 fixed = TRUE)
    expect_error(patchSize(cur_sce, convex = "test"),
                 regexp = "'convex' must be a single logical.",
                 fixed = TRUE)
})
