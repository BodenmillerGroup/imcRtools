test_that("testInteractions function works", {
    library(cytomapper)
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    ################################ classic ###################################
    expect_silent(cur_out <- testInteractions(pancreasSCE, 
                                                group_by = "ImageNb", 
                                                label = "CellType",
                                                colPairName = "knn_interaction_graph",
                                                iter = 100))
    
    expect_silent(cur_out <- testInteractions(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100,
                                              BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_silent(cur_out_2 <- testInteractions(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100,
                                                BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_identical(cur_out, cur_out_2)
    
    # Check numerical values
    expect_equal(cur_out$ct, c(1.8235294, 0.4705882, 0.7058824, 0.8750000, 0.2500000, 
                 1.8750000, 0.1250000, 0.1704545, 2.7045455, 2.3111111, 0.5111111, 
                 0.1777778, 0.7837838, 1.5675676, 0.6486486, 0.2000000, 0.6250000, 
                 2.1750000, NA, NA, NA, NA, 2.3636364, 0.6363636, NA, 0.7049180, 2.2950820), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_gt, c(0.00990099, 0.01980198, 1.00000000, 0.02970297, 0.53465347, 0.96039604, 
                                 1.00000000, 0.90099010, 0.00990099, 0.00990099, 1.00000000, 1.00000000, 
                                 0.99009901, 0.00990099, 1.00000000, 1.00000000, 0.99009901, 0.00990099,
                                 NA, NA, NA, NA, 0.00990099, 1.00000000, NA, 1.00000000, 0.00990099), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_lt, c(1.00000000, 0.99009901, 0.00990099, 1.00000000, 0.80198020, 
                                 0.06930693, 0.00990099, 0.19801980, 1.00000000, 1.00000000, 
                                 0.00990099, 0.00990099, 0.01980198, 1.00000000, 0.00990099, 
                                 0.00990099, 0.01980198, 1.00000000, NA, NA, NA, NA, 
                                 1.00000000, 0.00990099, NA, 0.00990099, 1.00000000), 
                 tolerance = 0.00001)
    
    # Check against countInteractions
    expect_silent(cur_sn <- countInteractions(pancreasSCE, 
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
    
    expect_silent(cur_out_2 <- testInteractions(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "classic",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100, p_threshold = 0.5,
                                                BPPARAM = SerialParam(RNGseed = 123)))    
    
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
    expect_silent(cur_out <- testInteractions(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "histocat",
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100,
                                              BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_silent(cur_out_2 <- testInteractions(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "histocat",
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100,
                                                BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_identical(cur_out, cur_out_2)
    
    # Check numerical values
    expect_equal(cur_out$ct, c(1.937500, 1.142857, 1.333333, 1.166667, 1.000000, 
                               1.875000, 1.222222, 1.363636, 2.735632, 2.476190, 
                               1.352941, 1.142857, 1.450000, 1.757576, 1.263158, 
                               1.333333, 1.666667, 2.636364, NA, NA, NA, NA,
                               2.363636, 1.312500, NA, 1.653846, 2.456140), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_gt, c(0.02970297, 0.19801980, 1.00000000, 0.43564356, 
                                 0.71287129, 0.98019802, 0.19801980, 0.00990099, 
                                 0.00990099, 0.00990099, 0.51485149, 1.00000000, 
                                 0.57425743, 0.00990099, 0.93069307, 0.88118812, 
                                 0.00990099, 0.00990099, NA, NA, NA, NA, 0.00990099,
                                 1.00000000, NA, 0.83168317, 0.00990099), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_lt, c(0.98019802, 0.81188119, 0.00990099, 0.57425743, 0.94059406, 
                                 0.02970297, 0.82178218, 1.00000000, 1.00000000, 
                                 1.00000000, 0.51485149, 0.00990099, 0.43564356, 
                                 1.00000000, 0.07920792, 0.15841584, 1.00000000, 
                                 1.00000000, NA, NA, NA, NA, 1.00000000, 0.00990099,
                                 NA, 0.18811881, 1.00000000), 
                 tolerance = 0.00001)
    
    # Check against countInteractions
    expect_silent(cur_sn <- countInteractions(pancreasSCE, 
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
    
    expect_silent(cur_out_2 <- testInteractions(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "histocat",
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100, p_threshold = 0.5,
                                                BPPARAM = SerialParam(RNGseed = 123)))    
    
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
    expect_silent(cur_out <- testInteractions(pancreasSCE, 
                                                  group_by = "ImageNb", 
                                                  label = "CellType",
                                                  method = "patch",
                                                  patch_size = 3,
                                                  colPairName = "knn_interaction_graph",
                                                  iter = 100,
                                              BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_silent(cur_out_2 <- testInteractions(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "patch",
                                                    patch_size = 3,
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100,
                                                BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_identical(cur_out, cur_out_2)
    
    # Check numerical values
    expect_equal(cur_out$ct, c(0.23529412, 0.00000000, 0.00000000, 0.00000000, 
                               0.00000000, 0.25000000, 0.00000000, 0.00000000, 
                               0.81818182, 0.55555556, 0.02222222, 0.00000000, 
                               0.02702703, 0.05405405, 0.00000000, 0.00000000, 
                               0.12500000, 0.52500000, NA, NA, NA, NA, 0.51515152, 
                               0.00000000, NA, 0.06557377, 0.57377049), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_gt, c(0.00990099, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 
                                 0.98019802, 1.00000000, 1.00000000, 0.00990099, 0.00990099, 
                                 0.63366337, 1.00000000, 0.82178218, 0.21782178, 1.00000000, 
                                 1.00000000, 0.00990099, 0.00990099, NA, NA, NA, NA, 0.00990099,
                                 1.00000000, NA, 0.99009901, 0.00990099), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_lt, c(1.00000000, 1.00000000, 0.00990099, 0.97029703, 1.00000000, 
                                 0.18811881, 0.80198020, 0.97029703, 1.00000000, 1.00000000,
                                 0.64356436, 0.22772277, 0.43564356, 0.94059406, 0.28712871, 
                                 0.14851485, 1.00000000, 1.00000000, NA, NA, NA, NA,
                                 1.00000000, 0.00990099, NA, 0.03960396, 1.00000000), 
                 tolerance = 0.00001)
    
    # Check against countInteractions
    expect_silent(cur_sn <- countInteractions(pancreasSCE, 
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
    
    expect_silent(cur_out_2 <- testInteractions(pancreasSCE, 
                                                    group_by = "ImageNb", 
                                                    label = "CellType",
                                                    method = "patch",
                                                    patch_size = 3,
                                                    colPairName = "knn_interaction_graph",
                                                    iter = 100, p_threshold = 0.5,
                                                BPPARAM = SerialParam(RNGseed = 123)))    
    
    expect_equal(cur_out_2$group_by, cur_sn$group_by)
    expect_equal(cur_out_2$from_label, cur_sn$from_label)
    expect_equal(cur_out_2$to_label, cur_sn$to_label)
    expect_equal(cur_out_2$ct, cur_sn$ct)
    
    cur_test <- cur_out_2[!is.na(cur_out_2$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.5, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)  
    
    # Corner case settings
    # one cell of a given cell-type and no neighbors
    cur_sce <- pancreasSCE
    cur_sce$CellType[123] <- "test"
    cur_sce <- buildSpatialGraph(cur_sce, img_id = "ImageNb", 
                                 type = "expansion", threshold = 7)
    
    plotSpatial(cur_sce, node_color_by = "CellType", img_id = "ImageNb", 
                draw_edges = TRUE, colPairName = "expansion_interaction_graph")
    
    expect_silent(cur_out <- testInteractions(cur_sce, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              colPairName = "expansion_interaction_graph",
                                              iter = 100,
                                              BPPARAM = SerialParam(RNGseed = 123)))
    
    # Check numerical values
    expect_equal(cur_out$ct, c(1.0000000, 0.0000000, 0.0000000, NA, 0.0000000, 0.0000000, 
                               0.0000000,        NA, 0.0000000, 0.0000000, 1.0588235,
                               NA,        NA,        NA,        NA,        NA, 
                               0.7500000, 0.3750000, 0.0000000, 0.0000000, 0.4285714, 0.2857143,
                               0.4285714, 0.0000000, 0.0000000, 0.4285714, 0.5714286, 0.0000000,
                               0.0000000, 0.0000000, 0.0000000, 0.0000000,        NA,
                               NA,        NA,        NA,        NA, 1.0526316, 0.2105263,        
                               NA,        NA, 0.1818182, 1.0000000,        NA,
                               NA,        NA,        NA,        NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_gt, c(0.00990099, 1.00000000, 1.00000000,         NA, 1.00000000, 
                                 1.00000000, 1.00000000,         NA, 1.00000000, 1.00000000,
                                 0.00990099,         NA,         NA,         NA,         NA,
                                 NA, 0.07920792, 0.41584158, 1.00000000, 1.00000000,
                                 0.42574257, 0.63366337, 0.43564356, 1.00000000, 1.00000000, 
                                 0.37623762, 0.18811881, 1.00000000, 1.00000000, 1.00000000,
                                 1.00000000, 1.00000000,        NA,         NA,         NA,
                                 NA,         NA, 0.00990099, 1.00000000,         NA,
                                 NA, 1.00000000, 0.00990099,         NA,         NA,         NA,         
                                 NA,         NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_lt, c(1, 0.782178217821782, 0.0297029702970297, NA, 0.782178217821782, 
                                 0.97029702970297, 0.247524752475248, NA, 0.0297029702970297, 
                                 0.247524752475248, 1, NA, NA, NA, NA, NA, 0.95049504950495, 0.623762376237624, 
                                 0.0594059405940594, 0.930693069306931, 0.613861386138614, 0.534653465346535, 
                                 0.613861386138614, 0.940594059405941, 0.0594059405940594, 0.683168316831683, 
                                 0.861386138613861, 0.95049504950495, 0.930693069306931, 0.940594059405941, 
                                 0.95049504950495, 1, NA, NA, NA, NA, NA, 1, 0.0099009900990099, 
                                 NA, NA, 0.0099009900990099, 1, NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    
    
    # Check against countInteractions
    expect_silent(cur_sn <- countInteractions(cur_sce, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              colPairName = "expansion_interaction_graph"))
    expect_equal(cur_out$group_by, cur_sn$group_by)
    expect_equal(cur_out$from_label, cur_sn$from_label)
    expect_equal(cur_out$to_label, cur_sn$to_label)
    expect_equal(cur_out$ct, cur_sn$ct)
    
    cur_test <- cur_out[!is.na(cur_out$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.01, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)
    
    ## histocat
    expect_silent(cur_out <- testInteractions(cur_sce, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              method = "histocat",
                                              colPairName = "expansion_interaction_graph",
                                              iter = 100,
                                              BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_equal(cur_out$ct, c(1.000000, 0.000000, 0.000000,       NA, 0.000000, 0.000000, 
                               0.000000,       NA, 0.000000, 0.000000, 1.058824,       NA,
                               NA,       NA,       NA,       NA, 1.000000, 1.000000, 
                               0.000000, 0.000000, 1.000000, 1.000000, 1.000000, 0.000000,
                               0.000000, 1.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 
                               0.000000,       NA,       NA,       NA,       NA,
                               NA, 1.250000, 1.000000,       NA,       NA, 1.000000, 
                               1.100000,       NA,       NA,       NA,       NA,       NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_gt, c(0.158415841584158, 1, 1, NA, 1, 1, 1, NA, 1, 1, 0.495049504950495, 
                                 NA, NA, NA, NA, NA, 0.792079207920792, 0.920792079207921, 1, 
                                 1, 0.920792079207921, 0.702970297029703, 0.910891089108911, 1, 
                                 1, 0.910891089108911, 0.742574257425743, 1, 1, 1, 1, 1, NA, NA, 
                                 NA, NA, NA, 0.099009900990099, 1, NA, NA, 1, 0.465346534653465, 
                                 NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_lt, c(0.99009900990099, 0.782178217821782, 0.0297029702970297, NA, 
                                 0.782178217821782, 0.97029702970297, 0.247524752475248, NA, 0.0297029702970297, 
                                 0.247524752475248, 0.574257425742574, NA, NA, NA, NA, NA, 0.900990099009901, 
                                 0.930693069306931, 0.0594059405940594, 0.930693069306931, 0.900990099009901, 
                                 0.99009900990099, 0.96039603960396, 0.940594059405941, 0.0594059405940594, 
                                 0.95049504950495, 0.910891089108911, 0.95049504950495, 0.930693069306931, 
                                 0.940594059405941, 0.95049504950495, 1, NA, NA, NA, NA, NA, 0.940594059405941, 
                                 0.297029702970297, NA, NA, 0.297029702970297, 0.544554455445545, 
                                 NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    
    # Check against countInteractions
    expect_silent(cur_sn <- countInteractions(cur_sce, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              method = "histocat",
                                              colPairName = "expansion_interaction_graph"))
    expect_equal(cur_out$group_by, cur_sn$group_by)
    expect_equal(cur_out$from_label, cur_sn$from_label)
    expect_equal(cur_out$to_label, cur_sn$to_label)
    expect_equal(cur_out$ct, cur_sn$ct)
    
    cur_test <- cur_out[!is.na(cur_out$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.01, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)
    
    ## patch
    expect_silent(cur_out <- testInteractions(cur_sce, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              method = "patch",
                                              patch_size = 1,
                                              colPairName = "expansion_interaction_graph",
                                              iter = 100,
                                              BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_equal(cur_out$ct, c(1.0000000, 0.0000000, 0.0000000,        NA, 0.0000000, 
                               0.0000000, 0.0000000,        NA, 0.0000000, 0.0000000, 1.0000000,
                               NA,        NA,        NA,        NA,        NA, 0.7500000, 
                               0.3750000, 0.0000000, 0.0000000, 0.4285714, 0.2857143,
                               0.4285714, 0.0000000, 0.0000000, 0.4285714, 0.5714286, 0.0000000, 
                               0.0000000, 0.0000000, 0.0000000, 0.0000000,        NA,
                               NA,        NA,        NA,        NA, 0.8421053, 0.2105263,        NA,
                               NA, 0.1818182, 0.9090909,        NA,
                               NA,        NA,        NA,        NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_gt, c(0.0099009900990099, 1, 1, NA, 1, 1, 1, NA, 1, 1, 0.0198019801980198, 
                                 NA, NA, NA, NA, NA, 0.0594059405940594, 0.396039603960396, 1, 
                                 1, 0.405940594059406, 0.633663366336634, 0.435643564356436, 1, 
                                 1, 0.376237623762376, 0.138613861386139, 1, 1, 1, 1, 1, NA, NA, 
                                 NA, NA, NA, 0.0198019801980198, 1, NA, NA, 1, 0.0099009900990099, 
                                 NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_lt, c(1, 0.782178217821782, 0.0297029702970297, NA, 0.782178217821782, 
                                 0.97029702970297, 0.247524752475248, NA, 0.0297029702970297, 
                                 0.247524752475248, 1, NA, NA, NA, NA, NA, 0.97029702970297, 0.653465346534653, 
                                 0.0594059405940594, 0.930693069306931, 0.653465346534653, 0.534653465346535, 
                                 0.613861386138614, 0.940594059405941, 0.0594059405940594, 0.683168316831683, 
                                 0.900990099009901, 0.95049504950495, 0.930693069306931, 0.940594059405941, 
                                 0.95049504950495, 1, NA, NA, NA, NA, NA, 0.99009900990099, 0.0099009900990099, 
                                 NA, NA, 0.0099009900990099, 1, NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    
    # Check against countInteractions
    expect_silent(cur_sn <- countInteractions(cur_sce, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              method = "patch",
                                              patch_size = 1,
                                              colPairName = "expansion_interaction_graph"))
    expect_equal(cur_out$group_by, cur_sn$group_by)
    expect_equal(cur_out$from_label, cur_sn$from_label)
    expect_equal(cur_out$to_label, cur_sn$to_label)
    expect_equal(cur_out$ct, cur_sn$ct)
    
    cur_test <- cur_out[!is.na(cur_out$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.01, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)
    
    # Fail
    expect_error(testInteractions("test"),
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "test", label = "CellType", 
                                      colPairName = "knn_interaction_graph"),
                 regexp = "'group_by' not in colData(object).",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = 1, label = "CellType", 
                                      colPairName = "knn_interaction_graph"),
                 regexp = "'group_by' must be a single string.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "test", 
                                      colPairName = "knn_interaction_graph"),
                 regexp = "'label' not in colData(object).",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = 1, 
                                      colPairName = "knn_interaction_graph"),
                 regexp = "'label' must be a single string.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "test"),
                 regexp = "'colPairName' not in colPairNames(object).",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", colPairName = 1),
                 regexp = "'colPairName' must be a single string.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                       method = "patch"),
                 regexp = "When method = 'patch', please specify 'patch_size'.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                       method = "patch", patch_size = "test"),
                 regexp = "'patch_size' must be a single numeric.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      iter = "test"),
                 regexp = "'iter' must be a single positive numeric.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      iter = c(1,2)),
                 regexp = "'iter' must be a single positive numeric.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      iter = -20),
                 regexp = "'iter' must be a single positive numeric.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      p_threshold = "test"),
                 regexp = "'p_threshold' must be a single numeric between 0 and 1.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      p_threshold = c(1,2)),
                 regexp = "'p_threshold' must be a single numeric between 0 and 1.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                      colPairName = "knn_interaction_graph",
                                      p_threshold = 3),
                 regexp = "'p_threshold' must be a single numeric between 0 and 1.",
                 fixed = TRUE)    
    

})
