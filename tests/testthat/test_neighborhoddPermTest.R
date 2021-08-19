test_that("neighborhoodPermTest function works", {
    library(cytomapper)
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    ################################ classic ###################################
    expect_silent(cur_out <- neighborhoodPermTest(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  colPairName = "knn_interaction_graph"))
    
    set.seed(123)
    expect_silent(cur_out <- neighborhoodPermTest(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100))
    
    set.seed(123)
    expect_silent(cur_out_2 <- neighborhoodPermTest(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100))
    
    expect_identical(cur_out, cur_out_2)
    
    # Check against summarizeNeigborhood
    expect_silent(cur_sn <- summarizeNeighborhood(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph"))
    expect_equal(cur_out$group_by, cur_sn$group_by)
    expect_equal(cur_out$from_label, cur_sn$from_label)
    expect_equal(cur_out$to_label, cur_sn$to_label)
    expect_equal(cur_out$ct, cur_sn$ct)
    
    cur_test <- cur_out[!is.na(cur_out$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.01, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)
    
    set.seed(123)
    expect_silent(cur_out_2 <- neighborhoodPermTest(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100, p_threshold = 0.5))    
    
    expect_equal(cur_out_2$group_by, cur_sn$group_by)
    expect_equal(cur_out_2$from_label, cur_sn$from_label)
    expect_equal(cur_out_2$to_label, cur_sn$to_label)
    expect_equal(cur_out_2$ct, cur_sn$ct)

    cur_test <- cur_out_2[!is.na(cur_out_2$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.5, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)  
    
    ################################ histocat ###################################
    set.seed(123)
    expect_silent(cur_out <- neighborhoodPermTest(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "histocat",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100))
    
    set.seed(123)
    expect_silent(cur_out_2 <- neighborhoodPermTest(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "histocat",
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100))
    
    expect_identical(cur_out, cur_out_2)
    
    # Check against summarizeNeigborhood
    expect_silent(cur_sn <- summarizeNeighborhood(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "histocat",
                                                  colPairName = "knn_interaction_graph"))
    expect_equal(cur_out$group_by, cur_sn$group_by)
    expect_equal(cur_out$from_label, cur_sn$from_label)
    expect_equal(cur_out$to_label, cur_sn$to_label)
    expect_equal(cur_out$ct, cur_sn$ct)
    
    cur_test <- cur_out[!is.na(cur_out$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.01, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)
    
    set.seed(123)
    expect_silent(cur_out_2 <- neighborhoodPermTest(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "histocat",
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100, p_threshold = 0.5))    
    
    expect_equal(cur_out_2$group_by, cur_sn$group_by)
    expect_equal(cur_out_2$from_label, cur_sn$from_label)
    expect_equal(cur_out_2$to_label, cur_sn$to_label)
    expect_equal(cur_out_2$ct, cur_sn$ct)
    
    cur_test <- cur_out_2[!is.na(cur_out_2$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.5, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)  
    
    
    ################################ patch ###################################
    set.seed(123)
    expect_silent(cur_out <- neighborhoodPermTest(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "patch",
                                                  patch_size = 3,
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100))
    
    set.seed(123)
    expect_silent(cur_out_2 <- neighborhoodPermTest(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "patch",
                                                    patch_size = 3,
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100))
    
    expect_identical(cur_out, cur_out_2)
    
    # Check against summarizeNeigborhood
    expect_silent(cur_sn <- summarizeNeighborhood(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "patch",
                                                  patch_size = 3,
                                                  colPairName = "knn_interaction_graph"))
    expect_equal(cur_out$group_by, cur_sn$group_by)
    expect_equal(cur_out$from_label, cur_sn$from_label)
    expect_equal(cur_out$to_label, cur_sn$to_label)
    expect_equal(cur_out$ct, cur_sn$ct)
    
    cur_test <- cur_out[!is.na(cur_out$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.01, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)
    
    set.seed(123)
    expect_silent(cur_out_2 <- neighborhoodPermTest(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "patch",
                                                    patch_size = 3,
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100, p_threshold = 0.5))    
    
    expect_equal(cur_out_2$group_by, cur_sn$group_by)
    expect_equal(cur_out_2$from_label, cur_sn$from_label)
    expect_equal(cur_out_2$to_label, cur_sn$to_label)
    expect_equal(cur_out_2$ct, cur_sn$ct)
    
    cur_test <- cur_out_2[!is.na(cur_out_2$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.5, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)  
    
    # Fail
    expect_error(neighborhoodPermTest("test"),
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "test", label = "CellType", 
                                      colPairName = "knn_interaction_graph"),
                 regexp = "'group_by' not in colData(object).",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = 1, label = "CellType", 
                                      colPairName = "knn_interaction_graph"),
                 regexp = "'group_by' must be a single string.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "test", 
                                      colPairName = "knn_interaction_graph"),
                 regexp = "'label' not in colData(object).",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = 1, 
                                      colPairName = "knn_interaction_graph"),
                 regexp = "'label' must be a single string.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "test"),
                 regexp = "'colPairName' not in colPairNames(object).",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", colPairName = 1),
                 regexp = "'colPairName' must be a single string.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                       method = "patch"),
                 regexp = "When method = 'patch', please specify 'patch_size'.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                       method = "patch", patch_size = "test"),
                 regexp = "'patch_size' must be a single numeric.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      iter = "test"),
                 regexp = "'iter' must be a single positive numeric.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      iter = c(1,2)),
                 regexp = "'iter' must be a single positive numeric.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      iter = -20),
                 regexp = "'iter' must be a single positive numeric.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      p_threshold = "test"),
                 regexp = "'p_threshold' must be a single numeric between 0 and 1.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      p_threshold = c(1,2)),
                 regexp = "'p_threshold' must be a single numeric between 0 and 1.",
                 fixed = TRUE)
    expect_error(neighborhoodPermTest(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      p_threshold = 3),
                 regexp = "'p_threshold' must be a single numeric between 0 and 1.",
                 fixed = TRUE)    
    

})
