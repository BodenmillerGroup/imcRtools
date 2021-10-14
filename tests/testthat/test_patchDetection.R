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
    
    expect_true(is(cur_sce, "SingleCellExperiment"))
    expect_true("patch_id" %in% names(colData(cur_sce)))
    expect_equal(unique(cur_sce$patch_id), c(NA, "1", "2", "3", "4", "5", "6", "7", "8"))
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType %in% c("celltype_B", "celltype_A"),
                                            colPairName = "expansion_interaction_graph"))
    expect_equal(unique(cur_sce$patch_id), c(NA, "1", "2", "3"))
    
    plotSpatial(cur_sce, img_id = "ImageNb", node_color_by = "patch_id")
    
    expect_silent(cur_sce <- patchDetection(pancreasSCE, 
                                            patch_cells = pancreasSCE$CellType == "celltype_B",
                                            colPairName = "expansion_interaction_graph",
                                            expand_by = 1))
    
    # Other graph constructors
    
    # Spatial Experiment
                     
})
