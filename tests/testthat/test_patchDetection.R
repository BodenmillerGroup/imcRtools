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
    
    
    # Spatial Experiment
    
    # Error
    
                     
})
