test_that("testInteractions function works", {
    library(cytomapper)
    library(BiocParallel)
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
    
    expect_equal(cur_out$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out$to_label, as.character(cur_sn$to_label))
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
    
    expect_equal(cur_out_2$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out_2$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out_2$to_label, as.character(cur_sn$to_label))
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
    expect_equal(cur_out$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out$to_label, as.character(cur_sn$to_label))
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
    
    expect_equal(cur_out_2$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out_2$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out_2$to_label, as.character(cur_sn$to_label))
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
    expect_equal(cur_out$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out$to_label, as.character(cur_sn$to_label))
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
    
    expect_equal(cur_out_2$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out_2$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out_2$to_label, as.character(cur_sn$to_label))
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
    expect_equal(cur_out$ct, c(0.235294117647059, 0, 0, NA, 0, 0, 0, NA, 0, 0, 0.204545454545455, 
                               NA, NA, NA, NA, NA, 0.133333333333333, 0.0666666666666667, 0, 
                               0, 0.0810810810810811, 0.0540540540540541, 0.0810810810810811, 
                               0, 0, 0.0769230769230769, 0.102564102564103, 0, 0, 0, 0, 0, NA, 
                               NA, NA, NA, NA, 0.303030303030303, 0.0606060606060606, NA, NA, 
                               0.0655737704918033, 0.360655737704918, NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_gt, c(0.0198019801980198, 1, 1, NA, 1, 1, 1, NA, 1, 1, 0.108910891089109, 
                                 NA, NA, NA, NA, NA, 0.227722772277228, 0.504950495049505, 1, 
                                 1, 0.504950495049505, 0.702970297029703, 0.455445544554455, 1, 
                                 1, 0.455445544554455, 0.376237623762376, 1, 1, 1, 1, 1, NA, NA, 
                                 NA, NA, NA, 0.108910891089109, 1, NA, NA, 1, 0.0198019801980198, 
                                 NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_lt, c(1, 0.782178217821782, 0.0297029702970297, NA, 0.782178217821782, 
                                 0.97029702970297, 0.247524752475248, NA, 0.0297029702970297, 
                                 0.247524752475248, 0.99009900990099, NA, NA, NA, NA, NA, 0.940594059405941, 
                                 0.732673267326733, 0.0594059405940594, 0.930693069306931, 0.732673267326733, 
                                 0.693069306930693, 0.762376237623762, 0.940594059405941, 0.0594059405940594, 
                                 0.762376237623762, 0.940594059405941, 0.95049504950495, 0.930693069306931, 
                                 0.940594059405941, 0.95049504950495, 1, NA, NA, NA, NA, NA, 0.96039603960396, 
                                 0.0099009900990099, NA, NA, 0.0099009900990099, 1, NA, NA, NA, 
                                 NA, NA), 
                 tolerance = 0.00001)
    
    
    # Check against countInteractions
    expect_silent(cur_sn <- countInteractions(cur_sce, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              colPairName = "expansion_interaction_graph"))
    expect_equal(cur_out$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out$to_label, as.character(cur_sn$to_label))
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
    expect_equal(cur_out$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out$to_label, as.character(cur_sn$to_label))
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
    
    expect_equal(cur_out$ct, c(0.235294117647059, 0, 0, NA, 0, 0, 0, NA, 0, 0, 0.193181818181818, 
                               NA, NA, NA, NA, NA, 0.133333333333333, 0.0666666666666667, 0, 
                               0, 0.0810810810810811, 0.0540540540540541, 0.0810810810810811, 
                               0, 0, 0.0769230769230769, 0.102564102564103, 0, 0, 0, 0, 0, NA, 
                               NA, NA, NA, NA, 0.242424242424242, 0.0606060606060606, NA, NA, 
                               0.0655737704918033, 0.327868852459016, NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_gt, c(0.0099009900990099, 1, 1, NA, 1, 1, 1, NA, 1, 1, 0.108910891089109, 
                                 NA, NA, NA, NA, NA, 0.198019801980198, 0.495049504950495, 1, 
                                 1, 0.495049504950495, 0.702970297029703, 0.455445544554455, 1, 
                                 1, 0.445544554455446, 0.326732673267327, 1, 1, 1, 1, 1, NA, NA, 
                                 NA, NA, NA, 0.168316831683168, 1, NA, NA, 1, 0.0099009900990099, 
                                 NA, NA, NA, NA, NA), 
                 tolerance = 0.00001)
    expect_equal(cur_out$p_lt, c(1, 0.782178217821782, 0.0297029702970297, NA, 0.782178217821782, 
                                 0.97029702970297, 0.247524752475248, NA, 0.0297029702970297, 
                                 0.247524752475248, 0.96039603960396, NA, NA, NA, NA, NA, 0.940594059405941, 
                                 0.752475247524752, 0.0594059405940594, 0.930693069306931, 0.792079207920792, 
                                 0.693069306930693, 0.762376237623762, 0.940594059405941, 0.0594059405940594, 
                                 0.772277227722772, 0.940594059405941, 0.95049504950495, 0.930693069306931, 
                                 0.940594059405941, 0.95049504950495, 1, NA, NA, NA, NA, NA, 0.900990099009901, 
                                 0.0099009900990099, NA, NA, 0.0099009900990099, 1, NA, NA, NA, 
                                 NA, NA), 
                 tolerance = 0.00001)
    
    # Check against countInteractions
    expect_silent(cur_sn <- countInteractions(cur_sce, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              method = "patch",
                                              patch_size = 1,
                                              colPairName = "expansion_interaction_graph"))
    expect_equal(cur_out$group_by, as.character(cur_sn$group_by))
    expect_equal(cur_out$from_label, as.character(cur_sn$from_label))
    expect_equal(cur_out$to_label, as.character(cur_sn$to_label))
    expect_equal(cur_out$ct, cur_sn$ct)
    
    cur_test <- cur_out[!is.na(cur_out$ct),]
    expect_equal(rowMin(as.matrix(cur_test[,c("p_gt", "p_lt")])), cur_test$p)
    expect_equal(cur_test$p_gt < cur_test$p_lt, cur_test$interaction)
    expect_equal(cur_test$p < 0.01, cur_test$sig)
    expect_equal(cur_test$sig * sign(cur_test$interaction - 0.5), cur_test$sigval)
    
    # Return samples
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "expansion",
                                     threshold = 20)
    expect_silent(cur_out <- testInteractions(pancreasSCE, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              method = "classic",
                                              iter = 100,
                                              colPairName = "expansion_interaction_graph",
                                              return_samples = TRUE,
                                              BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_silent(cur_out_2 <- testInteractions(pancreasSCE, 
                                              group_by = "ImageNb", 
                                              label = "CellType",
                                              method = "classic",
                                              iter = 100,
                                              colPairName = "expansion_interaction_graph",
                                              BPPARAM = SerialParam(RNGseed = 123)))
    
    expect_equal(cur_out$group_by, as.character(cur_out_2$group_by))
    expect_equal(cur_out$from_label, as.character(cur_out_2$from_label))
    expect_equal(cur_out$to_label, as.character(cur_out_2$to_label))
    expect_equal(cur_out$ct, cur_out_2$ct)
    expect_equal(cur_out$p_lt, cur_out_2$p_lt)
    expect_equal(cur_out$p_gt, cur_out_2$p_gt)
    expect_equal(cur_out$interaction, cur_out_2$interaction)
    expect_equal(cur_out$p, cur_out_2$p)
    expect_equal(cur_out$sig, cur_out_2$sig)
    expect_equal(cur_out$sigval, cur_out_2$sigval)
    
    expect_equal(dim(cur_out), c(27, 110))
    expect_equal(cur_out$iter_1, c(1.29411764705882, 0.352941176470588, 7.29411764705882, 0.75, 
                                   0.25, 7.125, 1.40909090909091, 0.647727272727273, 7.34090909090909, 
                                   3.64444444444444, 3.22222222222222, 2.95555555555556, 3.91891891891892, 
                                   3.08108108108108, 3.37837837837838, 3.325, 3.125, 3.3, 0, 0, 
                                   0, 0, 5.6969696969697, 5.51515151515152, 0, 5.9672131147541, 
                                   4.81967213114754))
    expect_equal(cur_out$iter_13, c(1.05882352941176, 0.470588235294118, 7.64705882352941, 1, 0.5, 
                                    6.5, 1.47727272727273, 0.590909090909091, 7.29545454545454, 4.31111111111111, 
                                    2.88888888888889, 2.75555555555556, 3.51351351351351, 3.67567567567568, 
                                    3.16216216216216, 3.1, 2.925, 3.6, 0, 0, 0, 0, 6.36363636363636, 
                                    5.03030303030303, 0, 5.44262295081967, 5.14754098360656))
    expect_equal(cur_out$iter_26, c(1.29411764705882, 0.882352941176471, 7.58823529411765, 1.875, 
                                    0.25, 7.125, 1.46590909090909, 0.647727272727273, 7.02272727272727, 
                                    3.91111111111111, 3.42222222222222, 3.11111111111111, 4.16216216216216, 
                                    2.54054054054054, 3.13513513513514, 3.5, 2.9, 3.15, 0, 0, 0, 
                                    0, 6.06060606060606, 5.09090909090909, 0, 5.50819672131148, 5.34426229508197
    ))
    expect_equal(cur_out$iter_67, c(1.64705882352941, 0.529411764705882, 6.29411764705882, 1.125, 
                                    0.25, 6.75, 1.21590909090909, 0.613636363636364, 7.65909090909091, 
                                    4.53333333333333, 2.86666666666667, 3.44444444444444, 3.48648648648649, 
                                    2.32432432432432, 3.13513513513514, 3.875, 2.9, 3.15, 0, 0, 0, 
                                    0, 5.96969696969697, 5.33333333333333, 0, 5.77049180327869, 4.91803278688525
    ))
    expect_equal(cur_out$iter_99, c(1.41176470588235, 0, 7.35294117647059, 0, 1.75, 8, 1.42045454545455, 
                                    0.727272727272727, 7.13636363636364, 3.91111111111111, 3.13333333333333, 
                                    3.53333333333333, 3.81081081081081, 2.86486486486486, 3.08108108108108, 
                                    3.975, 2.85, 2.65, 0, 0, 0, 0, 6.15151515151515, 5.31818181818182, 
                                    0, 5.75409836065574, 4.75409836065574))
    
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
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                  colPairName = "knn_interaction_graph",
                                  return_samples = 1),
                 regexp = "'return_samples' must be a single logical.",
                 fixed = TRUE)
    expect_error(testInteractions(pancreasSCE, group_by = "ImageNb", label = "CellType", 
                                  colPairName = "knn_interaction_graph",
                                  return_samples = c("test", "test2")),
                 regexp = "'return_samples' must be a single logical.",
                 fixed = TRUE)
    

})

test_that("testInteractions function works if cells are not grouped by image", {
    library(cytomapper)
    library(BiocParallel)
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "expansion",
                                     threshold = 10)
    
    plotSpatial(pancreasSCE, img_id = "ImageNb", draw_edges = TRUE, 
                colPairName = "expansion_interaction_graph", 
                node_color_by = "CellType", scales = "free")
    
    out <- testInteractions(pancreasSCE, group_by = "ImageNb", 
                            label = "CellType", method = "classic", 
                            iter = 200, colPairName = "expansion_interaction_graph", 
                            BPPARAM = SerialParam(RNGseed = 111))
    
    set.seed(123)
    cur_sce <- pancreasSCE[,sample(ncol(pancreasSCE))]
    
    plotSpatial(cur_sce, img_id = "ImageNb", draw_edges = TRUE, 
                colPairName = "expansion_interaction_graph", 
                node_color_by = "CellType", scales = "free")
    
    out2 <- testInteractions(cur_sce, group_by = "ImageNb", 
                            label = "CellType", method = "classic", 
                            iter = 200, 
                            colPairName = "expansion_interaction_graph",
                            BPPARAM = SerialParam(RNGseed = 111))
    
    ## Sampling is different between grouped and non-grouped cells
    ## Therefore the p values are close but not identical
    expect_equal(out$group_by, out2$group_by)
    expect_equal(out$from_label, out2$from_label)
    expect_equal(out$to_label, out2$to_label)
    expect_equal(out$ct, out2$ct)
    expect_equal(out$p_gt, out2$p_gt, tolerance = 0.1)
    expect_equal(out$p_lt, out2$p_lt, tolerance = 0.1)
    expect_equal(out$interaction, out2$interaction)
    expect_equal(out$p, out2$p, tolerance = 0.1)
    expect_equal(out$sig, out2$sig)
    expect_equal(out$sigval, out2$sigval)
    
    out <- testInteractions(pancreasSCE, group_by = "ImageNb", 
                            label = "CellType", method = "histocat", 
                            iter = 400, colPairName = "expansion_interaction_graph", 
                            BPPARAM = SerialParam(RNGseed = 111))
    
    out2 <- testInteractions(cur_sce, group_by = "ImageNb", 
                             label = "CellType", method = "histocat", 
                             iter = 400, 
                             colPairName = "expansion_interaction_graph",
                             BPPARAM = SerialParam(RNGseed = 111))
    
    expect_equal(out$group_by, out2$group_by)
    expect_equal(out$from_label, out2$from_label)
    expect_equal(out$to_label, out2$to_label)
    expect_equal(out$ct, out2$ct)
    expect_equal(out$p_gt, out2$p_gt, tolerance = 0.1)
    expect_equal(out$p_lt, out2$p_lt, tolerance = 0.1)
    expect_equal(out$interaction, out2$interaction)
    expect_equal(out$p, out2$p, tolerance = 0.1)
    #expect_equal(out$sig, out2$sig)
    #expect_equal(out$sigval, out2$sigval)
    
    out <- testInteractions(pancreasSCE, group_by = "ImageNb", 
                            label = "CellType", method = "patch", patch_size = 2,
                            iter = 300, colPairName = "expansion_interaction_graph", 
                            BPPARAM = SerialParam(RNGseed = 111))
    
    out2 <- testInteractions(cur_sce, group_by = "ImageNb", 
                             label = "CellType", method = "patch", patch_size = 2,
                             iter = 300, 
                             colPairName = "expansion_interaction_graph",
                             BPPARAM = SerialParam(RNGseed = 111))
    
    expect_equal(out$group_by, out2$group_by)
    expect_equal(out$from_label, out2$from_label)
    expect_equal(out$to_label, out2$to_label)
    expect_equal(out$ct, out2$ct)
    expect_equal(out$p_gt, out2$p_gt, tolerance = 0.1)
    expect_equal(out$p_lt, out2$p_lt, tolerance = 0.1)
    expect_equal(out$interaction, out2$interaction)
    expect_equal(out$p, out2$p, tolerance = 0.1)
    #expect_equal(out$sig, out2$sig)
    #expect_equal(out$sigval, out2$sigval)
    })
