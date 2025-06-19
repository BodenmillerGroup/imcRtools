test_that("distToCells works",{
  library(cytomapper)
  data("pancreasSCE")

  ################################ min ###################################
  # works when cell types present and with negative distances returned
  expect_message(cur_sce <- distToCells(object = pancreasSCE,
                                        x_cells = pancreasSCE$CellType == "celltype_B",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "min",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(class(cur_sce$distToCells) == "numeric")
  expect_true(min(cur_sce$distToCells) < 0)
  expect_true(sum(cur_sce$distToCells < 0) == sum(pancreasSCE$CellType == "celltype_B"))
  
  # enforces that a cells not of interest and a cell of interest do not have the same exact coordinates
  key_b <- paste(cur_sce[, cur_sce$CellType == "celltype_B"]$ImageName, cur_sce[, cur_sce$CellType == "celltype_B"]$Pos_X, cur_sce[, cur_sce$CellType == "celltype_B"]$Pos_Y, sep = "_")
  key_others <- paste(cur_sce[, cur_sce$CellType != "celltype_B"]$ImageName, cur_sce[, cur_sce$CellType != "celltype_B"]$Pos_X, cur_sce[, cur_sce$CellType != "celltype_B"]$Pos_Y, sep = "_")
  expect_true(length(intersect(key_b, key_others)) == 0)
  
  
  # Check numerical values
  expect_equal(unname(cur_sce$distToCells), c(79.788514, 20.857982, 53.368632, 62.713162, 19.560774, 41.448102, 62.462909, 68.450475, 31.699103, 54.291193, 10.058905, 46.570420, 10.015750,
                                      26.991639, 67.997428, 58.723036, 44.426290, 37.227309, 8.063107, 75.276820, 50.494166, -8.063107, 32.241412, 16.468327, 40.245508, 22.590794,
                                      62.698965, 8.458200, 74.717318, 52.895397, 8.385764, 23.644304, 43.704819, -8.385764, 30.637507, 37.041543, 68.431299, 58.295693, 14.364836,
                                      20.813568, 11.278759, 47.024721, 40.181783, 27.752307, 10.359226, 21.988126, 34.019184, 14.707679, 54.687570, 19.241309, 12.977059, 15.312316,
                                      42.197148, 8.949435, -8.949435, 66.681568, 27.796160, 18.739477, 22.029459, 37.334715, 22.518996, 52.515098, 11.601994, 66.141582, 16.125213,
                                      12.843852, 33.887107, 13.455102, 10.732585, 10.876951, 52.984788, 39.007635, 29.025695, 18.679265, -10.413550, 13.119897, 8.090589, 40.615390,
                                      47.265621, 23.254276, 13.716552, 57.794933, 26.927262, 17.825989, -8.090589, 37.227384, 8.909119, 53.916055, 15.965842, 18.554472, 22.880392,
                                      12.641206, 32.028055, -9.654366, 11.071430, 42.140034, 8.438556, 18.105918, -7.948255, 50.993988, 9.687375, 26.029224, 7.948255, 15.507842,
                                      14.496034, 16.470980, 32.213032, 9.190722, 44.430200, 22.904913, 9.436652, 37.564472, -9.190722, -12.910337, -10.600673, 10.513401, 15.610557,
                                      9.259873, -10.513401, 12.910337, 5.989048, 30.529686, 10.473304, 23.145159, -7.594326, -5.989048, 6.772751, 20.576832, -8.785245, -10.473304,
                                      -10.254327, 7.594326, -11.766886, -8.596174, 28.091328, 11.766886, -6.772751, -12.199488, 37.935954, -8.689778, 13.261366, 10.254327, 8.710711,
                                      -12.033538, 34.517442, 8.689778, -12.696385, 32.294570, 17.533721, -6.546034, -8.710711, -7.079958, 21.152374, 27.434023, 17.272200, 39.931658,
                                      23.705760, 6.546034, 6.861632, -10.260356, 14.630283, -8.866708, 26.171040, -6.861632, -13.335854, 33.213214, 22.572841, 25.816664, 10.531397,
                                      20.718661, 19.812147, 25.684989, -11.340031, 9.086057, 17.751145, 11.065690, 11.145864, 15.201760, 26.260594, 11.110730, 10.298755, 14.183898,
                                      -9.899502, 22.064641, 9.441764, -9.292279, -11.145864, 16.649749, 9.899502, 13.422591, 8.560659, -8.560659, 22.194993, 9.292279, 16.259552,
                                      17.388653, 18.876664, 8.719332, 11.814189, 19.147438, -10.300102, 11.414058, 7.699807, 15.151973, -7.699807, 14.767537, 17.681805, 8.860711,
                                      19.371295, 27.000206, 8.247212, -8.860711, 7.697301, 11.324157, 11.173417, -11.414058, 10.554889, 18.875406, -6.905692, -7.561787, 10.905143,
                                      -11.079098, 6.948386, 35.096189, 15.642309, 9.714982, 6.905692, 7.561787, 23.413358, 11.263664, 14.419926, -6.948386, -8.758330, 34.476882,
                                      -9.600399, 7.491662, -11.480456, -18.940514, -10.887572, -22.092775, -23.148566, -16.404099, -8.680828, -22.762279, -7.491662, -9.101670, -13.344017,
                                      -12.429175, -14.120853, 9.101670, -18.368374, 8.162161, -14.691913, -13.422799, 16.593679, -6.576539, -8.882700, -9.167578, -15.728945, -8.708470,
                                      6.576539, -11.189356, -20.543154, 7.616855, -18.566680, 8.708470, -13.159324, -23.342385, -7.616855, -12.370270, -15.567709, -6.933578, -7.852017,
                                      9.319316, 6.933578, -12.854097, 14.835364, 10.032582, -15.755407, -13.560830, -15.468886, -17.518757, -7.061241, -7.872385, 19.589301, 7.852017,
                                      8.030257, -9.400828, -6.685993, -8.588830, 15.229090, 12.196889, -23.506580, -14.180339, -11.017441, 8.588830, 16.856880, 24.159483, -19.166023,
                                      -6.531178, -9.428008, 6.685993, 18.854725, -8.827252, 6.531178, 19.869505, 8.827252, 25.368812, 28.841900, -21.495775, -9.266821, 9.064598,
                                      17.816143, 35.230210, 29.085772, -18.158671, 25.665566, 36.054209, -7.867960, 7.867960, 35.011391, 7.970550, -10.121218, -10.783356, 10.763131,
                                      42.537944, -14.673660, -11.795126, -7.970550, -9.963250, 21.561468, 42.260122, 10.872469, -7.003065, 29.467776, 10.121218, -13.828144, 7.003065,
                                      14.871748, -7.179175, -13.455011, 48.715264, 11.213500, 37.520319, 7.179175, 26.735469, -11.213500, -12.026209, 21.238797, 35.155338, 46.384431,
                                      15.379208, 24.752166, 9.317498, 19.900208, -11.012489, 23.229380, -9.911600, 11.012489, 29.874798, 13.596656, 29.695571), 
               tolerance = 0.00001)

  # works on cell types when present and no negative distances returned
  expect_message(cur_sce_2 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_B",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "min",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_2, "SingleCellExperiment"))
  expect_s4_class(cur_sce_2 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_2)))
  expect_true(class(cur_sce_2$distToCells) == "numeric")
  expect_true(min(cur_sce_2$distToCells) == 0)
  expect_true(sum(cur_sce_2$distToCells == 0) >= sum(pancreasSCE$CellType == "celltype_B"))

  expect_equal(cur_sce[,cur_sce$distToCells > 0]$distToCells,cur_sce_2[,cur_sce_2$distToCells > 0]$distToCells)

  expect_equal(length(cur_sce[,cur_sce$distToCells < 0]),length(cur_sce_2[,cur_sce_2$distToCells == 0]))

  # works on cell types when not present in some image and with negative distances returned
  expect_message(cur_sce_3 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_A",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "min",
                                          img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_3, "SingleCellExperiment"))
  expect_s4_class(cur_sce_3 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_3)))
  expect_true(is(cur_sce_3$distToCells, "numeric"))
  expect_true(sum(cur_sce_3$distToCells[cur_sce_3$ImageName != "J02_imc.tiff"] < 0) == sum(pancreasSCE$CellType == "celltype_A"))

  expect_true(any(is.na(cur_sce_3$distToCells)))
  expect_true(all(is.na(cur_sce_3[,cur_sce_3$ImageName == "J02_imc.tiff"]$distToCells)))
  expect_true(min(cur_sce_3[,!is.na(cur_sce_3$distToCells)]$distToCells)<0)

  expect_equal(length(cur_sce_3[,cur_sce_3$ImageName == "J02_imc.tiff"]$CellNb),sum(is.na(cur_sce_3$distToCells)))

  # works on cell types when not present in some images and no negative distances returned
  expect_message(cur_sce_4 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_A",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "min",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_4, "SingleCellExperiment"))
  expect_s4_class(cur_sce_4 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_4)))
  expect_true(is(cur_sce_4$distToCells, "numeric"))
  expect_true(sum(cur_sce_4$distToCells[cur_sce_4$ImageName != "J02_imc.tiff"] < 0) == 0)

  expect_true(any(is.na(cur_sce_4$distToCells)))
  expect_true(all(is.na(cur_sce_4[,cur_sce_4$ImageName == "J02_imc.tiff"]$distToCells)))
  expect_true(min(cur_sce_4[,!is.na(cur_sce_4$distToCells)]$distToCells) == 0)

  expect_equal(length(cur_sce_4[,cur_sce_4$ImageName == "J02_imc.tiff"]$CellNb),sum(is.na(cur_sce_4$distToCells)))

  # Spatial Experiment
  cur_spe <- SpatialExperiment:::.sce_to_spe(pancreasSCE, sample_id = as.character(pancreasSCE$ImageNb))
  spatialCoords(cur_spe) <- as.matrix(colData(pancreasSCE)[,c("Pos_X", "Pos_Y")])
  colData(cur_spe)[c("Pos_X", "Pos_Y")] <- NULL

  cur_spe_1 <- distToCells(cur_spe,
                           x_cells = cur_spe$CellType == "celltype_B",
                           coords = c("Pos_X","Pos_Y"),
                           metric = "min",
                           img_id = "ImageName")

  expect_true(is(cur_spe_1, "SingleCellExperiment"))
  expect_s4_class(cur_spe_1 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_spe_1)))
  expect_true(class(cur_spe_1$distToCells) == "numeric")
  expect_true(min(cur_spe_1$distToCells) < 0)
  
  # Check numerical values
  expect_equal(unname(cur_spe_1$distToCells),
               c(79.788514, 20.857982, 53.368632, 62.713162, 19.560774, 41.448102, 62.462909, 68.450475, 
                 31.699103, 54.291193, 10.058905, 46.570420, 10.015750, 26.991639, 67.997428, 58.723036, 
                 44.426290, 37.227309, 8.063107, 75.276820, 50.494166, -8.063107, 32.241412, 16.468327, 
                 40.245508, 22.590794, 62.698965, 8.458200, 74.717318, 52.895397, 8.385764, 23.644304, 
                 43.704819, -8.385764, 30.637507, 37.041543, 68.431299, 58.295693, 14.364836, 20.813568, 
                 11.278759, 47.024721, 40.181783, 27.752307, 10.359226, 21.988126, 34.019184, 14.707679, 
                 54.687570, 19.241309, 12.977059, 15.312316, 42.197148, 8.949435, -8.949435, 66.681568, 
                 27.796160, 18.739477, 22.029459, 37.334715, 22.518996, 52.515098, 11.601994, 66.141582, 
                 16.125213, 12.843852, 33.887107, 13.455102, 10.732585, 10.876951, 52.984788, 39.007635, 
                 29.025695, 18.679265, -10.413550, 13.119897, 8.090589, 40.615390, 47.265621, 23.254276, 
                 13.716552, 57.794933, 26.927262, 17.825989, -8.090589, 37.227384, 8.909119, 53.916055, 
                 15.965842, 18.554472, 22.880392, 12.641206, 32.028055, -9.654366, 11.071430, 42.140034, 
                 8.438556, 18.105918, -7.948255, 50.993988, 9.687375, 26.029224, 7.948255, 15.507842, 
                 14.496034, 16.470980, 32.213032, 9.190722, 44.430200, 22.904913, 9.436652, 37.564472, 
                 -9.190722, -12.910337, -10.600673, 10.513401, 15.610557, 9.259873, -10.513401, 12.910337, 
                 5.989048, 30.529686, 10.473304, 23.145159, -7.594326, -5.989048, 6.772751, 20.576832, 
                 -8.785245, -10.473304, -10.254327, 7.594326, -11.766886, -8.596174, 28.091328, 11.766886, 
                 -6.772751, -12.199488, 37.935954, -8.689778, 13.261366, 10.254327, 8.710711, -12.033538, 
                 34.517442, 8.689778, -12.696385, 32.294570, 17.533721, -6.546034, -8.710711, -7.079958, 
                 21.152374, 27.434023, 17.272200, 39.931658, 23.705760, 6.546034, 6.861632, -10.260356, 
                 14.630283, -8.866708, 26.171040, -6.861632, -13.335854, 33.213214, 22.572841, 25.816664, 
                 10.531397, 20.718661, 19.812147, 25.684989, -11.340031, 9.086057, 17.751145, 11.065690, 
                 11.145864, 15.201760, 26.260594, 11.110730, 10.298755, 14.183898, -9.899502, 22.064641, 
                 9.441764, -9.292279, -11.145864, 16.649749, 9.899502, 13.422591, 8.560659, -8.560659, 
                 22.194993, 9.292279, 16.259552, 17.388653, 18.876664, 8.719332, 11.814189, 19.147438, 
                 -10.300102, 11.414058, 7.699807, 15.151973, -7.699807, 14.767537, 17.681805, 8.860711, 
                 19.371295, 27.000206, 8.247212, -8.860711, 7.697301, 11.324157, 11.173417, -11.414058, 
                 10.554889, 18.875406, -6.905692, -7.561787, 10.905143, -11.079098, 6.948386, 35.096189, 
                 15.642309, 9.714982, 6.905692, 7.561787, 23.413358, 11.263664, 14.419926, -6.948386, -8.758330, 
                 34.476882, -9.600399, 7.491662, -11.480456, -18.940514, -10.887572, -22.092775, -23.148566, 
                 -16.404099, -8.680828, -22.762279, -7.491662, -9.101670, -13.344017, -12.429175, -14.120853, 
                 9.101670, -18.368374, 8.162161, -14.691913, -13.422799, 16.593679, -6.576539, -8.882700, 
                 -9.167578, -15.728945, -8.708470, 6.576539, -11.189356, -20.543154, 7.616855, -18.566680, 
                 8.708470, -13.159324, -23.342385, -7.616855, -12.370270, -15.567709, -6.933578, -7.852017, 
                 9.319316, 6.933578, -12.854097, 14.835364, 10.032582, -15.755407, -13.560830, -15.468886, 
                 -17.518757, -7.061241, -7.872385, 19.589301, 7.852017, 8.030257, -9.400828, -6.685993, 
                 -8.588830, 15.229090, 12.196889, -23.506580, -14.180339, -11.017441, 8.588830, 16.856880, 
                 24.159483, -19.166023, -6.531178, -9.428008, 6.685993, 18.854725, -8.827252, 6.531178, 19.869505, 
                 8.827252, 25.368812, 28.841900, -21.495775, -9.266821, 9.064598, 17.816143, 35.230210, 29.085772, 
                 -18.158671, 25.665566, 36.054209, -7.867960, 7.867960, 35.011391, 7.970550, -10.121218, -10.783356, 
                 10.763131, 42.537944, -14.673660, -11.795126, -7.970550, -9.963250, 21.561468, 42.260122, 10.872469, 
                 -7.003065, 29.467776, 10.121218, -13.828144, 7.003065, 14.871748, -7.179175, -13.455011, 48.715264, 
                 11.213500, 37.520319, 7.179175, 26.735469, -11.213500, -12.026209, 21.238797, 35.155338, 46.384431, 
                 15.379208, 24.752166, 9.317498, 19.900208, -11.012489, 23.229380, -9.911600, 11.012489, 29.874798, 
                 13.596656, 29.695571), 
               tolerance = 0.00001)

  # works on cell types when present and no negative distances returned
  expect_message(cur_spe_2 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_B",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "min",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_spe_2, "SingleCellExperiment"))
  expect_s4_class(cur_spe_2 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_spe_2)))
  expect_true(class(cur_spe_2$distToCells) == "numeric")
  expect_true(min(cur_spe_2$distToCells) == 0)

  expect_equal(cur_spe_1[,cur_spe_1$distToCells > 0]$distToCells,cur_spe_2[,cur_spe_2$distToCells > 0]$distToCells)

  expect_equal(length(cur_spe_1[,cur_spe_1$distToCells < 0]),length(cur_spe_2[,cur_spe_2$distToCells == 0]))

  # compare results from SingleCellExperiment and SpatialExperiment
  expect_equal(cur_sce$distToCells,cur_spe_1$distToCells)

  expect_equal(cur_sce_2$distToCells,cur_spe_2$distToCells)

  # Works when all cells of an image belong to one batch
  expect_message(cur_sce <- distToCells(object = pancreasSCE,
                                        x_cells = pancreasSCE$ImageName == "J02_imc.tiff",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "min",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(all(is.na(cur_sce$distToCells)))

  cur_sce$CellType[cur_sce$ImageName == "J02_imc.tiff"] <- "celltype_A"
  expect_message(cur_sce <- distToCells(object = cur_sce,
                                        x_cells = cur_sce$CellType == "celltype_A",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "min",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")
  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(all(is.na(cur_sce$distToCells[cur_sce$ImageName == "J02_imc.tiff"])))
  expect_true(all(!is.na(cur_sce$distToCells[cur_sce$ImageName != "J02_imc.tiff"])))

  ################################ max ###################################
  # works when cell types present and with negative distances returned
  expect_message(cur_sce <- distToCells(object = pancreasSCE,
                                        x_cells = pancreasSCE$CellType == "celltype_B",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "max",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(class(cur_sce$distToCells) == "numeric")
  expect_true(min(cur_sce$distToCells) < 0)
  expect_true(sum(cur_sce$distToCells < 0) == sum(pancreasSCE$CellType == "celltype_B"))
  
  # enforces that a cells not of interest and a cell of interest do not have the same exact coordinates
  key_b <- paste(cur_sce[, cur_sce$CellType == "celltype_B"]$ImageName, cur_sce[, cur_sce$CellType == "celltype_B"]$Pos_X, cur_sce[, cur_sce$CellType == "celltype_B"]$Pos_Y, sep = "_")
  key_others <- paste(cur_sce[, cur_sce$CellType != "celltype_B"]$ImageName, cur_sce[, cur_sce$CellType != "celltype_B"]$Pos_X, cur_sce[, cur_sce$CellType != "celltype_B"]$Pos_Y, sep = "_")
  expect_true(length(intersect(key_b, key_others)) == 0)
  
  # Check numerical values
  expect_equal(unname(cur_sce$distToCells), 
               c(136.65424, 109.45536, 114.55055, 122.57908, 103.20744, 110.10950, 121.08675, 125.85645, 104.02530, 113.88156, 
                 99.48648, 107.32173, 96.27177, 95.51628, 122.83985, 115.17086, 104.39610, 97.98693, 90.77567, 128.31391, 107.80927, 
                 -125.40640, 93.06610, 87.52394, 98.85802, 84.94491, 114.72467, 80.30537, 123.87844, 105.84826, 83.25590, 84.48384, 
                 99.03982, -112.70969, 89.11108, 91.75206, 114.01386, 106.89030, 74.79366, 80.56622, 70.67018, 97.41946, 92.76868, 
                 82.70275, 72.13011, 77.56976, 85.08788, 62.23649, 99.13506, 76.10229, 69.66951, 66.48275, 88.16187, 60.35536, 
                 -83.70912, 108.34969, 75.48661, 68.67006, 63.55929, 79.51520, 67.52256, 93.37157, 50.06821, 111.98730, 59.86538, 
                 53.13621, 73.92562, 55.07534, 53.93810, 52.10275, 100.84196, 83.79838, 71.04568, 59.82518, -102.41230, 59.61156, 
                 53.92648, 94.85896, 104.09008, 68.19041, 59.95850, 114.37683, 87.97371, 78.19425, -112.02505, 102.05000, 64.05678, 
                 118.96672, 88.32270, 71.85520, 97.64075, 73.16095, 105.24512, -130.14058, 80.40444, 113.34043, 91.35054, 81.78349, 
                 -108.01429, 125.40640, 79.76830, 107.43902, 89.74841, 87.06216, 101.36620, 85.25133, 114.84130, 84.32018, 124.76685, 
                 110.11188, 98.95015, 120.65315, -136.65424, -117.19962, -124.59190, 116.01144, 117.81951, 109.50685, -112.61780, 123.38543, 
                 105.02567, 133.29062, 110.30807, 129.21410, -114.28904, -115.08516, 109.21176, 120.50851, -118.93194, -121.31788, -100.79870, 
                 95.96305, -99.23385, -101.34053, 121.55165, 110.09809, -126.86419, -110.12156, 127.67964, -99.12318, 101.11575, 91.64309, 
                 96.47528, -113.70254, 118.97897, 83.16991, -99.84497, 111.82954, 93.07938, -94.79728, -121.06791, -110.43911, 85.57133, 
                 101.64536, 79.24731, 120.14710, 92.83363, 74.51657, 76.69385, -111.74206, 71.39877, -89.43099, 84.15705, -96.76329, -116.63364, 
                 115.26425, 90.48126, 107.58291, 64.36820, 96.01323, 69.94381, 75.83017, -85.37164, 76.76023, 98.96759, 80.55351, 87.66133, 
                 80.53120, 110.93248, 88.27672, 65.70795, 94.20606, -83.77524, 67.76377, 94.28646, -91.65969, -102.83113, 72.65571, 75.70841, 
                 98.06972, 81.04490, -100.84991, 107.26495, 83.38916, 78.93739, 105.20016, 81.16837, 95.51904, 90.96825, 86.96895, -110.43348, 
                 109.75437, 88.69014, 102.92720, -106.29110, 87.03683, 94.47550, 101.08877, 95.03576, 115.90965, 100.55947, -118.24822, 94.02984, 
                 103.44229, 101.79673, -133.29062, 97.80124, 111.19450, -120.36616, -110.78990, 109.73760, -130.14253, 104.46286, 126.86419, 
                 109.24756, 106.27451, 117.19962, 107.47124, 114.94604, 113.54808, 113.97164, -113.67717, -117.81951, 126.76258, -127.92093, 
                 105.15307, -128.14760, -113.72931, -121.56970, -129.20328, -140.03783, -125.35490, -115.25028, -134.21105, -115.33109, 
                 -119.74376, -119.31641, -109.70900, -128.54707, 94.47143, -104.10310, 95.21246, -133.98505, -102.94967, 98.56684, -103.41747, 
                 -109.03012, -128.38267, -111.57787, -119.24319, 86.43528, -109.31940, -97.01776, 95.46974, -101.98786, 103.98350, -88.91075, 
                 -94.89587, -111.01974, -107.99267, -103.48775, -87.50252, -108.82267, 113.41930, 73.83941, -111.48394, 106.31868, 99.73997, 
                 -97.36220, -80.40135, -86.13588, -103.53626, -83.56312, -91.56783, 110.72034, 88.52006, 76.30119, -93.77473, -86.04469, -79.89423, 
                 82.11458, 91.13685, -112.86446, -75.80021, -74.22879, 70.42659, 94.93553, 104.90421, -102.27092, -76.76421, -92.14173, 86.42590, 
                 86.95476, -67.48553, 84.65533, 81.50350, 68.78545, 96.32017, 102.79873, -104.40624, -85.38980, 78.44790, 78.76383, 107.58587, 
                 92.37917, -95.92974, 88.02283, 102.82469, -79.85229, 86.79890, 98.00726, 83.83443, -103.51643, -113.41930, 82.38212, 111.82921, 
                 -91.44376, -85.47864, -81.42985, -81.68047, 95.91262, 111.27608, 90.53437, -97.79017, 105.07726, 122.06475, -89.34922, 117.16862, 
                 124.32938, -96.05414, -89.32267, 123.24252, 101.97731, 116.60195, 117.11714, 111.46082, -95.27058, -96.05173, 131.74314, 119.56759, 
                 128.14760, 125.59791, 114.46212, 116.46141, 113.37011, -101.69945, 133.34071, -105.15307, 108.95504, 140.03783, 127.42944, 123.81012), 
               tolerance = 0.00001)

  # works on cell types when present and no negative distances returned
  expect_message(cur_sce_2 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_B",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "max",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_2, "SingleCellExperiment"))
  expect_s4_class(cur_sce_2 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_2)))
  expect_true(class(cur_sce_2$distToCells) == "numeric")

  expect_equal(cur_sce[,!pancreasSCE$CellType == "celltype_B"]$distToCells,cur_sce_2[,!pancreasSCE$CellType == "celltype_B"]$distToCells)

  expect_equal(sum(cur_sce$distToCells < 0),sum(pancreasSCE$CellType == "celltype_B"))

  # works on cell types when not present in some image and with negative distances returned
  expect_message(cur_sce_3 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_A",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "max",
                                          img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_3, "SingleCellExperiment"))
  expect_s4_class(cur_sce_3 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_3)))
  expect_true(is(cur_sce_3$distToCells, "numeric"))
  expect_true(sum(cur_sce_3[,cur_sce_3$ImageName != "J02_imc.tiff"]$distToCells < 0) == sum(pancreasSCE$CellType == "celltype_A"))

  expect_true(any(is.na(cur_sce_3$distToCells)))
  expect_true(all(is.na(cur_sce_3[,cur_sce_3$ImageName == "J02_imc.tiff"]$distToCells)))
  expect_true(min(cur_sce_3[,!is.na(cur_sce_3$distToCells)]$distToCells)<0)

  expect_equal(length(cur_sce_3[,cur_sce_3$ImageName == "J02_imc.tiff"]$CellNb),sum(is.na(cur_sce_3$distToCells)))

  # works on cell types when not present in some images and no negative distances returned
  expect_message(cur_sce_4 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_A",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "max",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_4, "SingleCellExperiment"))
  expect_s4_class(cur_sce_4 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_4)))
  expect_true(is(cur_sce_4$distToCells, "numeric"))
  expect_true(sum(cur_sce_4[,cur_sce_4$ImageName != "J02_imc.tiff"]$distToCells < 0) == 0)

  expect_true(any(is.na(cur_sce_4$distToCells)))
  expect_true(all(is.na(cur_sce_4[,cur_sce_4$ImageName == "J02_imc.tiff"]$distToCells)))
  expect_true(min(cur_sce_4[,!is.na(cur_sce_4$distToCells)]$distToCells) > 0)

  expect_equal(length(cur_sce_4[,cur_sce_4$ImageName == "J02_imc.tiff"]$CellNb),sum(is.na(cur_sce_4$distToCells)))

  # Spatial Experiment
  cur_spe <- SpatialExperiment:::.sce_to_spe(pancreasSCE, sample_id = as.character(pancreasSCE$ImageNb))
  spatialCoords(cur_spe) <- as.matrix(colData(pancreasSCE)[,c("Pos_X", "Pos_Y")])
  colData(cur_spe)[c("Pos_X", "Pos_Y")] <- NULL

  cur_spe_1 <- distToCells(cur_spe,
                           x_cells = cur_spe$CellType == "celltype_B",
                           coords = c("Pos_X","Pos_Y"),
                           metric = "max",
                           img_id = "ImageName")

  expect_true(is(cur_spe_1, "SingleCellExperiment"))
  expect_s4_class(cur_spe_1 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_spe_1)))
  expect_true(class(cur_spe_1$distToCells) == "numeric")
  expect_true(min(cur_spe_1$distToCells) < 0)
  
  # Check numerical values
  expect_equal(unname(cur_spe_1$distToCells), 
               c(136.65424, 109.45536, 114.55055, 122.57908, 103.20744, 110.10950, 121.08675, 125.85645, 104.02530, 113.88156, 
                 99.48648, 107.32173, 96.27177, 95.51628, 122.83985, 115.17086, 104.39610, 97.98693, 90.77567, 128.31391, 107.80927, 
                 -125.40640, 93.06610, 87.52394, 98.85802, 84.94491, 114.72467, 80.30537, 123.87844, 105.84826, 83.25590, 84.48384, 
                 99.03982, -112.70969, 89.11108, 91.75206, 114.01386, 106.89030, 74.79366, 80.56622, 70.67018, 97.41946, 92.76868, 
                 82.70275, 72.13011, 77.56976, 85.08788, 62.23649, 99.13506, 76.10229, 69.66951, 66.48275, 88.16187, 60.35536, 
                 -83.70912, 108.34969, 75.48661, 68.67006, 63.55929, 79.51520, 67.52256, 93.37157, 50.06821, 111.98730, 59.86538, 
                 53.13621, 73.92562, 55.07534, 53.93810, 52.10275, 100.84196, 83.79838, 71.04568, 59.82518, -102.41230, 59.61156, 
                 53.92648, 94.85896, 104.09008, 68.19041, 59.95850, 114.37683, 87.97371, 78.19425, -112.02505, 102.05000, 64.05678, 
                 118.96672, 88.32270, 71.85520, 97.64075, 73.16095, 105.24512, -130.14058, 80.40444, 113.34043, 91.35054, 81.78349, 
                 -108.01429, 125.40640, 79.76830, 107.43902, 89.74841, 87.06216, 101.36620, 85.25133, 114.84130, 84.32018, 124.76685, 
                 110.11188, 98.95015, 120.65315, -136.65424, -117.19962, -124.59190, 116.01144, 117.81951, 109.50685, -112.61780, 
                 123.38543, 105.02567, 133.29062, 110.30807, 129.21410, -114.28904, -115.08516, 109.21176, 120.50851, -118.93194, 
                 -121.31788, -100.79870, 95.96305, -99.23385, -101.34053, 121.55165, 110.09809, -126.86419, -110.12156, 127.67964, 
                 -99.12318, 101.11575, 91.64309, 96.47528, -113.70254, 118.97897, 83.16991, -99.84497, 111.82954, 93.07938, -94.79728, 
                 -121.06791, -110.43911, 85.57133, 101.64536, 79.24731, 120.14710, 92.83363, 74.51657, 76.69385, -111.74206, 71.39877, 
                 -89.43099, 84.15705, -96.76329, -116.63364, 115.26425, 90.48126, 107.58291, 64.36820, 96.01323, 69.94381, 75.83017, 
                 -85.37164, 76.76023, 98.96759, 80.55351, 87.66133, 80.53120, 110.93248, 88.27672, 65.70795, 94.20606, -83.77524, 67.76377, 
                 94.28646, -91.65969, -102.83113, 72.65571, 75.70841, 98.06972, 81.04490, -100.84991, 107.26495, 83.38916, 78.93739, 
                 105.20016, 81.16837, 95.51904, 90.96825, 86.96895, -110.43348, 109.75437, 88.69014, 102.92720, -106.29110, 87.03683, 
                 94.47550, 101.08877, 95.03576, 115.90965, 100.55947, -118.24822, 94.02984, 103.44229, 101.79673, -133.29062, 97.80124, 
                 111.19450, -120.36616, -110.78990, 109.73760, -130.14253, 104.46286, 126.86419, 109.24756, 106.27451, 117.19962, 107.47124, 
                 114.94604, 113.54808, 113.97164, -113.67717, -117.81951, 126.76258, -127.92093, 105.15307, -128.14760, -113.72931, 
                 -121.56970, -129.20328, -140.03783, -125.35490, -115.25028, -134.21105, -115.33109, -119.74376, -119.31641, -109.70900, 
                 -128.54707, 94.47143, -104.10310, 95.21246, -133.98505, -102.94967, 98.56684, -103.41747, -109.03012, -128.38267, -111.57787, 
                 -119.24319, 86.43528, -109.31940, -97.01776, 95.46974, -101.98786, 103.98350, -88.91075, -94.89587, -111.01974, -107.99267, 
                 -103.48775, -87.50252, -108.82267, 113.41930, 73.83941, -111.48394, 106.31868, 99.73997, -97.36220, -80.40135, -86.13588, 
                 -103.53626, -83.56312, -91.56783, 110.72034, 88.52006, 76.30119, -93.77473, -86.04469, -79.89423, 82.11458, 91.13685, 
                 -112.86446, -75.80021, -74.22879, 70.42659, 94.93553, 104.90421, -102.27092, -76.76421, -92.14173, 86.42590, 86.95476, 
                 -67.48553, 84.65533, 81.50350, 68.78545, 96.32017, 102.79873, -104.40624, -85.38980, 78.44790, 78.76383, 107.58587, 
                 92.37917, -95.92974, 88.02283, 102.82469, -79.85229, 86.79890, 98.00726, 83.83443, -103.51643, -113.41930, 82.38212, 
                 111.82921, -91.44376, -85.47864, -81.42985, -81.68047, 95.91262, 111.27608, 90.53437, -97.79017, 105.07726, 122.06475, 
                 -89.34922, 117.16862, 124.32938, -96.05414, -89.32267, 123.24252, 101.97731, 116.60195, 117.11714, 111.46082, -95.27058, 
                 -96.05173, 131.74314, 119.56759, 128.14760, 125.59791, 114.46212, 116.46141, 113.37011, -101.69945, 133.34071, -105.15307, 
                 108.95504, 140.03783, 127.42944, 123.81012), 
               tolerance = 0.00001)

  # works on cell types when present and no negative distances returned
  expect_message(cur_spe_2 <- distToCells(object = cur_spe,
                                          x_cells = cur_spe$CellType == "celltype_B",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "max",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_spe_2, "SingleCellExperiment"))
  expect_s4_class(cur_spe_2 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_spe_2)))
  expect_true(class(cur_spe_2$distToCells) == "numeric")
  expect_true(min(cur_spe_2$distToCells) > 0)

  expect_equal(cur_spe_1[,!cur_spe$CellType == "celltype_B"]$distToCells,cur_spe_2[,!cur_spe$CellType == "celltype_B"]$distToCells)

  expect_equal(sum(cur_spe_1$distToCells < 0),sum(cur_spe$CellType == "celltype_B"))

  # compare results from SingleCellExperiment and SpatialExperiment
  expect_equal(cur_sce$distToCells,cur_spe_1$distToCells)

  expect_equal(cur_sce_2$distToCells,cur_spe_2$distToCells)

  # Works when all cells of an image belong to one batch
  expect_message(cur_sce <- distToCells(object = pancreasSCE,
                                        x_cells = pancreasSCE$ImageName == "J02_imc.tiff",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "max",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(all(is.na(cur_sce$distToCells)))

  cur_sce$CellType[cur_sce$ImageName == "J02_imc.tiff"] <- "celltype_A"
  expect_message(cur_sce <- distToCells(object = cur_sce,
                                        x_cells = cur_sce$CellType == "celltype_A",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "max",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")
  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(all(is.na(cur_sce$distToCells[cur_sce$ImageName == "J02_imc.tiff"])))
  expect_true(all(!is.na(cur_sce$distToCells[cur_sce$ImageName != "J02_imc.tiff"])))


  ################################ mean ###################################
  # works when cell types present and with negative distances returned
  expect_message(cur_sce <- distToCells(object = pancreasSCE,
                                        x_cells = pancreasSCE$CellType == "celltype_B",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "mean",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(class(cur_sce$distToCells) == "numeric")
  expect_true(min(cur_sce$distToCells) < 0)
  expect_true(sum(cur_sce$distToCells < 0) == sum(pancreasSCE$CellType == "celltype_B"))
  
  # enforces that a cells not of interest and a cell of interest do not have the same exact coordinates
  key_b <- paste(cur_sce[, cur_sce$CellType == "celltype_B"]$ImageName, cur_sce[, cur_sce$CellType == "celltype_B"]$Pos_X, cur_sce[, cur_sce$CellType == "celltype_B"]$Pos_Y, sep = "_")
  key_others <- paste(cur_sce[, cur_sce$CellType != "celltype_B"]$ImageName, cur_sce[, cur_sce$CellType != "celltype_B"]$Pos_X, cur_sce[, cur_sce$CellType != "celltype_B"]$Pos_Y, sep = "_")
  expect_true(length(intersect(key_b, key_others)) == 0)
  
  # Check numerical values
  expect_equal(unname(cur_sce$distToCells), 
               c(105.85771, 72.29898, 79.61456, 89.17125, 65.68919, 73.80932, 89.27000, 95.01068, 67.33008, 81.38977, 62.37970, 
                 73.24674, 59.27292, 59.60926, 94.34333, 85.74286, 72.14951, 64.84985, 54.24479, 101.00038, 78.10317, -67.60954, 
                 59.63898, 51.93169, 68.70304, 51.17622, 88.72235, 46.05889, 99.52372, 79.65491, 46.71290, 52.66989, 71.59700, 
                 -58.88041, 60.17709, 65.38024, 92.36616, 83.57423, 44.60594, 51.74329, 40.25658, 73.44848, 67.71925, 57.11611, 
                 39.97433, 52.25172, 61.60701, 36.07447, 78.73353, 42.81982, 44.97609, 37.24455, 67.40372, 35.33337, -45.05630, 
                 89.95882, 54.68452, 47.42234, 36.00308, 61.30392, 48.69543, 76.40332, 32.09296, 95.16136, 42.13856, 32.14515, 
                 56.18728, 32.53570, 36.38179, 32.35866, 82.10385, 64.22791, 49.64015, 35.43279, -53.69635, 33.68054, 30.97870, 
                 72.04224, 80.85516, 39.94478, 37.86204, 91.97714, 60.50739, 47.76106, -59.44954, 74.42890, 33.10569, 92.18078, 
                 55.14932, 35.93583, 64.54653, 34.78900, 73.44214, -73.59345, 42.62600, 82.97328, 53.64393, 40.85210, -57.26668, 
                 93.53869, 39.90361, 71.29469, 47.29299, 44.35473, 62.08366, 41.91028, 78.02857, 44.21130, 89.35978, 70.42524, 
                 55.11644, 83.52813, -78.06966, -69.38935, -71.66376, 60.00098, 60.92113, 52.38330, -63.22311, 72.86176, 49.39583, 
                 84.74269, 53.61752, 79.77662, -63.12164, -63.22975, 54.16982, 72.11223, -65.88714, -68.00709, -55.33135, 43.77503, 
                 -55.76341, -54.49485, 75.40090, 62.84994, -73.17916, -59.50082, 82.75139, -52.26550, 57.27340, 47.56781, 47.34267, 
                 -63.07772, 76.25835, 41.41397, -53.11539, 70.95936, 52.55127, -49.78442, -70.44314, -61.44641, 47.86127, 62.48563, 
                 43.71106, 80.88443, 54.82339, 38.42133, 38.05371, -63.55159, 41.20550, -47.73303, 51.08923, -52.80394, -69.05632, 
                 78.92698, 56.73416, 73.19959, 39.29958, 62.32742, 43.52542, 47.53616, -46.89710, 38.99505, 67.28870, 41.17965, 
                 44.11888, 53.77657, 79.74015, 60.02963, 41.53005, 48.37929, -45.21830, 48.24092, 68.34432, -49.00068, -56.43345, 
                 46.31064, 43.65253, 49.68298, 56.79609, -43.78897, 80.42736, 43.61117, 52.52128, 76.98450, 48.25630, 67.41313, 
                 46.09382, 48.30664, -59.04680, 56.04556, 59.03220, 73.81331, -46.01597, 55.03709, 49.72119, 51.19446, 53.76951, 
                 85.16472, 69.87880, -63.68038, 59.43716, 54.29873, 69.04244, -77.49842, 59.04159, 79.32374, -64.60122, -51.85864, 
                 58.89070, -73.27876, 63.42316, 94.10801, 75.27418, 60.92222, 62.52133, 70.70554, 81.90207, 63.56832, 76.80485, 
                 -59.72884, -59.93289, 93.07046, -70.80516, 53.98722, -86.48224, -74.76748, -70.62280, -72.26216, -77.93084, -70.93176, 
                 -70.05600, -74.90978, -76.17393, -79.25530, -65.16082, -71.88829, -68.27229, 51.43734, -68.09353, 45.88415, -70.69189, 
                 -64.34451, 58.03419, -60.82901, -58.19483, -64.76527, -57.09399, -60.39487, 41.91318, -70.36574, -62.17182, 51.35777, 
                 -64.85942, 58.39617, -56.39332, -59.81838, -71.67019, -52.34546, -65.84488, -52.22994, -50.62153, 68.35546, 38.14989, 
                 -72.79719, 64.13270, 59.13449, -61.42579, -50.38716, -53.74205, -66.48498, -47.09967, -46.74877, 69.57877, 51.29759, 
                 43.82828, -59.20948, -53.79296, -43.02041, 49.37090, 55.12922, -75.31273, -47.60338, -44.96983, 43.11979, 59.45852, 
                 68.07657, -66.79105, -48.31046, -59.11257, 38.35676, 55.04872, -43.58501, 38.30565, 50.29805, 41.10757, 62.70091, 
                 68.71254, -67.48719, -53.15903, 41.88122, 46.69859, 71.86027, 58.13642, -59.32100, 53.42071, 66.48652, -47.67620, 
                 43.35110, 61.42004, 45.21701, -62.57813, -70.98051, 48.42981, 74.30457, -53.37670, -49.39061, -43.31792, -46.60550, 
                 58.53749, 72.63991, 53.90519, -55.93295, 66.35028, 61.39281, -50.53594, 58.79055, 64.19635, -54.04325, -48.92265, 
                 83.02016, 62.82217, 76.37337, 60.61733, 71.28331, -46.22300, -49.04979, 70.92167, 78.59864, 86.74930, 67.73478, 
                 73.67375, 63.46815, 72.59989, -53.23699, 73.74178, -57.46298, 69.15286, 79.83211, 72.91781, 81.98807), 
               tolerance = 0.00001)

  # works on cell types when present and no negative distances returned
  expect_message(cur_sce_2 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_B",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "mean",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_2, "SingleCellExperiment"))
  expect_s4_class(cur_sce_2 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_2)))
  expect_true(class(cur_sce_2$distToCells) == "numeric")

  expect_equal(cur_sce[,!pancreasSCE$CellType == "celltype_B"]$distToCells,cur_sce_2[,!pancreasSCE$CellType == "celltype_B"]$distToCells)

  expect_equal(sum(cur_sce$distToCells < 0),sum(pancreasSCE$CellType == "celltype_B"))

  # works on cell types when not present in some image and with negative distances returned
  expect_message(cur_sce_3 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_A",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "mean",
                                          img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_3, "SingleCellExperiment"))
  expect_s4_class(cur_sce_3 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_3)))
  expect_true(is(cur_sce_3$distToCells, "numeric"))
  expect_true(sum(cur_sce_3[,cur_sce_3$ImageName != "J02_imc.tiff"]$distToCells < 0) == sum(pancreasSCE$CellType == "celltype_A"))

  expect_true(any(is.na(cur_sce_3$distToCells)))
  expect_true(all(is.na(cur_sce_3[,cur_sce_3$ImageName == "J02_imc.tiff"]$distToCells)))
  expect_true(min(cur_sce_3[,!is.na(cur_sce_3$distToCells)]$distToCells)<0)

  expect_equal(length(cur_sce_3[,cur_sce_3$ImageName == "J02_imc.tiff"]$CellNb),sum(is.na(cur_sce_3$distToCells)))

  # works on cell types when not present in some images and no negative distances returned
  expect_message(cur_sce_4 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_A",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "mean",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_4, "SingleCellExperiment"))
  expect_s4_class(cur_sce_4 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_4)))
  expect_true(is(cur_sce_4$distToCells, "numeric"))
  expect_true(sum(cur_sce_4[,cur_sce_4$ImageName != "J02_imc.tiff"]$distToCells < 0) == 0)

  expect_true(any(is.na(cur_sce_4$distToCells)))
  expect_true(all(is.na(cur_sce_4[,cur_sce_4$ImageName == "J02_imc.tiff"]$distToCells)))
  expect_true(min(cur_sce_4[,!is.na(cur_sce_4$distToCells)]$distToCells) > 0)

  expect_equal(length(cur_sce_4[,cur_sce_4$ImageName == "J02_imc.tiff"]$CellNb),sum(is.na(cur_sce_4$distToCells)))

  # Spatial Experiment
  cur_spe <- SpatialExperiment:::.sce_to_spe(pancreasSCE, sample_id = as.character(pancreasSCE$ImageNb))
  spatialCoords(cur_spe) <- as.matrix(colData(pancreasSCE)[,c("Pos_X", "Pos_Y")])
  colData(cur_spe)[c("Pos_X", "Pos_Y")] <- NULL

  cur_spe_1 <- distToCells(cur_spe,
                           x_cells = cur_spe$CellType == "celltype_B",
                           coords = c("Pos_X","Pos_Y"),
                           metric = "mean",
                           img_id = "ImageName")

  expect_true(is(cur_spe_1, "SingleCellExperiment"))
  expect_s4_class(cur_spe_1 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_spe_1)))
  expect_true(class(cur_spe_1$distToCells) == "numeric")
  expect_true(min(cur_spe_1$distToCells) < 0)
  
  # Check numerical values
  expect_equal(unname(cur_spe_1$distToCells),
               c(105.85771, 72.29898, 79.61456, 89.17125, 65.68919, 73.80932, 89.27000, 95.01068, 67.33008, 81.38977, 62.37970, 73.24674, 59.27292, 59.60926, 94.34333,
                 85.74286, 72.14951, 64.84985, 54.24479, 101.00038, 78.10317, -67.60954, 59.63898, 51.93169, 68.70304, 51.17622, 88.72235, 46.05889, 99.52372, 79.65491,
                 46.71290, 52.66989, 71.59700, -58.88041, 60.17709, 65.38024, 92.36616, 83.57423, 44.60594, 51.74329, 40.25658, 73.44848, 67.71925, 57.11611, 39.97433,
                 52.25172, 61.60701, 36.07447, 78.73353, 42.81982, 44.97609, 37.24455, 67.40372, 35.33337, -45.05630, 89.95882, 54.68452, 47.42234, 36.00308, 61.30392,
                 48.69543, 76.40332, 32.09296, 95.16136, 42.13856, 32.14515, 56.18728, 32.53570, 36.38179, 32.35866, 82.10385, 64.22791, 49.64015, 35.43279, -53.69635,
                 33.68054, 30.97870, 72.04224, 80.85516, 39.94478, 37.86204, 91.97714, 60.50739, 47.76106, -59.44954, 74.42890, 33.10569, 92.18078, 55.14932, 35.93583,
                 64.54653, 34.78900, 73.44214, -73.59345, 42.62600, 82.97328, 53.64393, 40.85210, -57.26668, 93.53869, 39.90361, 71.29469, 47.29299, 44.35473, 62.08366,
                 41.91028, 78.02857, 44.21130, 89.35978, 70.42524, 55.11644, 83.52813, -78.06966, -69.38935, -71.66376, 60.00098, 60.92113, 52.38330, -63.22311, 72.86176,
                 49.39583, 84.74269, 53.61752, 79.77662, -63.12164, -63.22975, 54.16982, 72.11223, -65.88714, -68.00709, -55.33135, 43.77503, -55.76341, -54.49485, 75.40090,
                 62.84994, -73.17916, -59.50082, 82.75139, -52.26550, 57.27340, 47.56781, 47.34267, -63.07772, 76.25835, 41.41397, -53.11539, 70.95936, 52.55127, -49.78442,
                 -70.44314, -61.44641, 47.86127, 62.48563, 43.71106, 80.88443, 54.82339, 38.42133, 38.05371, -63.55159, 41.20550, -47.73303, 51.08923, -52.80394, -69.05632,
                 78.92698, 56.73416, 73.19959, 39.29958, 62.32742, 43.52542, 47.53616, -46.89710, 38.99505, 67.28870, 41.17965, 44.11888, 53.77657, 79.74015, 60.02963,
                 41.53005, 48.37929, -45.21830, 48.24092, 68.34432, -49.00068, -56.43345, 46.31064, 43.65253, 49.68298, 56.79609, -43.78897, 80.42736, 43.61117, 52.52128,
                 76.98450, 48.25630, 67.41313, 46.09382, 48.30664, -59.04680, 56.04556, 59.03220, 73.81331, -46.01597, 55.03709, 49.72119, 51.19446, 53.76951, 85.16472,
                 69.87880, -63.68038, 59.43716, 54.29873, 69.04244, -77.49842, 59.04159, 79.32374, -64.60122, -51.85864, 58.89070, -73.27876, 63.42316, 94.10801, 75.27418,
                 60.92222, 62.52133, 70.70554, 81.90207, 63.56832, 76.80485, -59.72884, -59.93289, 93.07046, -70.80516, 53.98722, -86.48224, -74.76748, -70.62280, -72.26216,
                 -77.93084, -70.93176, -70.05600, -74.90978, -76.17393, -79.25530, -65.16082, -71.88829, -68.27229, 51.43734, -68.09353, 45.88415, -70.69189, -64.34451, 58.03419,
                 -60.82901, -58.19483, -64.76527, -57.09399, -60.39487, 41.91318, -70.36574, -62.17182, 51.35777, -64.85942, 58.39617, -56.39332, -59.81838, -71.67019, -52.34546,
                 -65.84488, -52.22994, -50.62153, 68.35546, 38.14989, -72.79719, 64.13270, 59.13449, -61.42579, -50.38716, -53.74205, -66.48498, -47.09967, -46.74877, 69.57877,
                 51.29759, 43.82828, -59.20948, -53.79296, -43.02041, 49.37090, 55.12922, -75.31273, -47.60338, -44.96983, 43.11979, 59.45852, 68.07657, -66.79105, -48.31046,
                 -59.11257, 38.35676, 55.04872, -43.58501, 38.30565, 50.29805, 41.10757, 62.70091, 68.71254, -67.48719, -53.15903, 41.88122, 46.69859, 71.86027, 58.13642,
                 -59.32100, 53.42071, 66.48652, -47.67620, 43.35110, 61.42004, 45.21701, -62.57813, -70.98051, 48.42981, 74.30457, -53.37670, -49.39061, -43.31792, -46.60550,
                 58.53749, 72.63991, 53.90519, -55.93295, 66.35028, 61.39281, -50.53594, 58.79055, 64.19635, -54.04325, -48.92265, 83.02016, 62.82217, 76.37337, 60.61733,
                 71.28331, -46.22300, -49.04979, 70.92167, 78.59864, 86.74930, 67.73478, 73.67375, 63.46815, 72.59989, -53.23699, 73.74178, -57.46298, 69.15286, 79.83211,
                 72.91781, 81.98807), 
               tolerance = 0.00001)

  # works on cell types when present and no negative distances returned
  expect_message(cur_spe_2 <- distToCells(object = cur_spe,
                                          x_cells = cur_spe$CellType == "celltype_B",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "mean",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_spe_2, "SingleCellExperiment"))
  expect_s4_class(cur_spe_2 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_spe_2)))
  expect_true(class(cur_spe_2$distToCells) == "numeric")
  expect_true(min(cur_spe_2$distToCells) > 0)

  expect_equal(cur_spe_1[,!cur_spe$CellType == "celltype_B"]$distToCells,cur_spe_2[,!cur_spe$CellType == "celltype_B"]$distToCells)

  expect_equal(sum(cur_spe_1$distToCells < 0),sum(cur_spe$CellType == "celltype_B"))

  # compare results from SingleCellExperiment and SpatialExperiment
  expect_equal(cur_sce$distToCells,cur_spe_1$distToCells)

  expect_equal(cur_sce_2$distToCells,cur_spe_2$distToCells)

  # Works when all cells of an image belong to one batch
  expect_message(cur_sce <- distToCells(object = pancreasSCE,
                                        x_cells = pancreasSCE$ImageName == "J02_imc.tiff",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "mean",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(all(is.na(cur_sce$distToCells)))

  cur_sce$CellType[cur_sce$ImageName == "J02_imc.tiff"] <- "celltype_A"
  expect_message(cur_sce <- distToCells(object = cur_sce,
                                        x_cells = cur_sce$CellType == "celltype_A",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "mean",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")
  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(all(is.na(cur_sce$distToCells[cur_sce$ImageName == "J02_imc.tiff"])))
  expect_true(all(!is.na(cur_sce$distToCells[cur_sce$ImageName != "J02_imc.tiff"])))


  ################################ median ###################################
  # works when cell types present and with negative distances returned
  expect_message(cur_sce <- distToCells(object = pancreasSCE,
                                        x_cells = pancreasSCE$CellType == "celltype_B",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "median",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(class(cur_sce$distToCells) == "numeric")
  expect_true(min(cur_sce$distToCells) < 0)
  expect_true(sum(cur_sce$distToCells < 0) == sum(pancreasSCE$CellType == "celltype_B"))
  
  # enforces that a cells not of interest and a cell of interest do not have the same exact coordinates
  key_b <- paste(cur_sce[, cur_sce$CellType == "celltype_B"]$ImageName, cur_sce[, cur_sce$CellType == "celltype_B"]$Pos_X, cur_sce[, cur_sce$CellType == "celltype_B"]$Pos_Y, sep = "_")
  key_others <- paste(cur_sce[, cur_sce$CellType != "celltype_B"]$ImageName, cur_sce[, cur_sce$CellType != "celltype_B"]$Pos_X, cur_sce[, cur_sce$CellType != "celltype_B"]$Pos_Y, sep = "_")
  expect_true(length(intersect(key_b, key_others)) == 0)
  
  # Check numerical values
  expect_equal(unname(cur_sce$distToCells),
               c(103.79708, 77.19487, 81.74023, 90.91419, 70.00210, 76.60087, 89.34814, 93.48272, 70.34545, 82.43675, 67.03657, 74.67073, 63.67169, 61.83572, 90.00359,
                 82.94776, 72.62156, 65.34178, 58.06859, 95.41740, 76.17895, -68.19041, 59.98736, 54.03205, 68.04943, 51.31451, 84.47693, 47.12413, 96.63302, 75.31743,
                 47.52159, 51.50512, 68.10810, -60.37985, 58.29240, 63.16882, 92.24195, 82.34067, 41.53038, 49.83299, 37.17257, 71.62839, 65.04276, 57.01429, 35.28390,
                 52.83621, 59.54155, 30.91311, 78.95270, 37.94634, 46.58579, 30.50557, 67.07761, 33.64414, -44.87616, 91.79820, 53.63460, 45.51329, 32.09087, 61.96482,
                 48.08096, 78.62163, 35.08889, 98.68567, 41.07923, 30.60649, 57.69076, 29.36967, 38.73110, 37.84127, 85.76051, 67.08190, 51.47189, 35.98709, -54.98000,
                 34.60094, 33.31931, 76.29741, 84.79908, 38.14804, 39.55889, 96.16675, 62.34194, 47.04409, -62.09536, 76.14043, 33.02666, 94.00564, 53.68852, 31.73132,
                 63.22703, 26.67201, 72.73931, -77.93963, 39.11328, 82.91384, 50.36465, 31.60422, -56.01568, 92.26648, 38.68745, 67.77602, 41.94134, 36.17816, 57.93474,
                 29.68316, 73.77379, 44.08081, 85.78152, 65.89571, 51.10218, 78.84303, -80.56622, -69.28504, -71.90990, 57.84604, 51.47220, 42.17414, -62.23792, 72.98064,
                 40.98654, 83.77524, 42.83898, 79.68969, -63.46598, -65.09480, 45.18274, 71.86100, -68.20659, -69.78288, -53.65702, 35.70765, -56.32207, -53.92707, 75.51740,
                 61.24951, -76.91300, -61.20788, 83.65608, -52.39311, 56.52574, 42.38686, 44.15033, -65.45875, 76.32325, 34.84570, -53.88628, 70.13359, 50.31694, -48.90561,
                 -71.56121, -62.40373, 44.82710, 61.05060, 40.77808, 81.06731, 52.35672, 33.41591, 34.32866, -64.06619, 39.38950, -48.12136, 51.74043, -52.23378, -69.81853,
                 78.23132, 53.79857, 72.25755, 40.73339, 60.65922, 43.76330, 47.55786, -46.33360, 38.63039, 66.36261, 39.85121, 41.35095, 54.76625, 80.35763, 61.44830,
                 43.58614, 45.35274, -43.61361, 50.88073, 72.04698, -49.60022, -56.37487, 48.96206, 42.82675, 48.90561, 59.53004, -39.93166, 85.53775, 41.55653, 56.46341,
                 79.78497, 50.17230, 67.60904, 44.77366, 47.45894, -62.74491, 57.73427, 61.68409, 74.56752, -44.63916, 58.97305, 49.70444, 51.26106, 57.01699, 86.40865,
                 73.61289, -65.20768, 63.42481, 55.80264, 74.97227, -80.45364, 64.28639, 83.00309, -67.30757, -50.91488, 60.46133, -76.38780, 70.86505, 97.19144, 81.65832,
                 67.42034, 66.74744, 77.86342, 88.41234, 68.08210, 84.76199, -57.80302, -59.69768, 98.95125, -71.39877, 53.28076, -90.53726, -72.39165, -68.50196, -72.27422,
                 -80.52008, -68.54955, -66.72530, -77.07195, -80.17271, -83.93850, -63.25588, -76.28748, -70.02824, 51.10953, -71.14010, 40.46037, -72.17862, -62.57266, 57.67305,
                 -55.90097, -55.39778, -66.67122, -56.42756, -61.20298, 39.14313, -74.94019, -65.25053, 50.15070, -67.58988, 63.21401, -57.96498, -62.43940, -75.79812, -51.70963,
                 -68.06333, -51.82955, -48.40215, 74.88232, 37.92051, -73.54655, 68.43147, 62.15394, -63.81999, -51.80399, -56.94498, -68.26764, -48.44406, -42.89581, 72.97069,
                 52.43054, 44.65180, -60.19297, -55.56035, -44.19917, 49.11250, 56.88895, -79.00633, -50.70163, -49.26376, 44.74917, 60.85039, 70.45090, -68.78874, -49.40788,
                 -59.14884, 32.55601, 54.20772, -47.36248, 33.61630, 50.39182, 41.55033, 62.29585, 67.79165, -68.78558, -53.99483, 42.54096, 46.58633, 70.58747, 58.49726,
                 -63.52971, 54.37175, 65.23877, -49.07492, 41.10413, 61.03113, 44.97764, -66.87476, -77.32278, 49.24589, 72.17169, -57.06608, -50.47158, -43.81847, -47.49533,
                 60.17817, 71.54758, 56.43004, -60.83717, 68.67465, 56.98453, -51.91766, 56.47203, 61.09518, -57.10227, -47.87032, 82.64731, 66.26489, 77.77807, 59.42314,
                 73.24015, -45.39070, -46.03222, 67.89339, 80.60991, 87.66147, 66.61122, 77.62832, 64.50831, 76.77484, -50.84583, 72.57423, -54.62405, 74.23006, 78.38152,
                 73.46824, 86.50457), 
               tolerance = 0.00001)

  # works on cell types when present and no negative distances returned
  expect_message(cur_sce_2 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_B",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "median",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_2, "SingleCellExperiment"))
  expect_s4_class(cur_sce_2 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_2)))
  expect_true(class(cur_sce_2$distToCells) == "numeric")

  expect_equal(cur_sce[,!pancreasSCE$CellType == "celltype_B"]$distToCells,cur_sce_2[,!pancreasSCE$CellType == "celltype_B"]$distToCells)

  expect_equal(sum(cur_sce$distToCells < 0),sum(pancreasSCE$CellType == "celltype_B"))

  # works on cell types when not present in some image and with negative distances returned
  expect_message(cur_sce_3 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_A",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "median",
                                          img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_3, "SingleCellExperiment"))
  expect_s4_class(cur_sce_3 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_3)))
  expect_true(is(cur_sce_3$distToCells, "numeric"))
  expect_true(sum(cur_sce_3[,cur_sce_3$ImageName != "J02_imc.tiff"]$distToCells < 0) == sum(pancreasSCE$CellType == "celltype_A"))

  expect_true(any(is.na(cur_sce_3$distToCells)))
  expect_true(all(is.na(cur_sce_3[,cur_sce_3$ImageName == "J02_imc.tiff"]$distToCells)))
  expect_true(min(cur_sce_3[,!is.na(cur_sce_3$distToCells)]$distToCells)<0)

  expect_equal(length(cur_sce_3[,cur_sce_3$ImageName == "J02_imc.tiff"]$CellNb),sum(is.na(cur_sce_3$distToCells)))

  # works on cell types when not present in some images and no negative distances returned
  expect_message(cur_sce_4 <- distToCells(object = pancreasSCE,
                                          x_cells = pancreasSCE$CellType == "celltype_A",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "median",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_sce_4, "SingleCellExperiment"))
  expect_s4_class(cur_sce_4 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce_4)))
  expect_true(is(cur_sce_4$distToCells, "numeric"))
  expect_true(sum(cur_sce_4[,cur_sce_4$ImageName != "J02_imc.tiff"]$distToCells < 0) == 0)

  expect_true(any(is.na(cur_sce_4$distToCells)))
  expect_true(all(is.na(cur_sce_4[,cur_sce_4$ImageName == "J02_imc.tiff"]$distToCells)))
  expect_true(min(cur_sce_4[,!is.na(cur_sce_4$distToCells)]$distToCells) > 0)

  expect_equal(length(cur_sce_4[,cur_sce_4$ImageName == "J02_imc.tiff"]$CellNb),sum(is.na(cur_sce_4$distToCells)))

  # Spatial Experiment
  cur_spe <- SpatialExperiment:::.sce_to_spe(pancreasSCE, sample_id = as.character(pancreasSCE$ImageNb))
  spatialCoords(cur_spe) <- as.matrix(colData(pancreasSCE)[,c("Pos_X", "Pos_Y")])
  colData(cur_spe)[c("Pos_X", "Pos_Y")] <- NULL

  cur_spe_1 <- distToCells(cur_spe,
                           x_cells = cur_spe$CellType == "celltype_B",
                           coords = c("Pos_X","Pos_Y"),
                           metric = "median",
                           img_id = "ImageName")

  expect_true(is(cur_spe_1, "SingleCellExperiment"))
  expect_s4_class(cur_spe_1 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_spe_1)))
  expect_true(class(cur_spe_1$distToCells) == "numeric")
  expect_true(min(cur_spe_1$distToCells) < 0)
  
  # Check numerical values
  expect_equal(unname(cur_spe_1$distToCells), 
               c(103.79708, 77.19487, 81.74023, 90.91419, 70.00210, 76.60087, 89.34814, 93.48272, 70.34545, 82.43675, 67.03657, 74.67073, 63.67169, 61.83572, 90.00359, 82.94776, 72.62156, 65.34178, 58.06859, 95.41740, 76.17895, -68.19041, 59.98736, 54.03205, 68.04943, 51.31451, 84.47693, 47.12413, 96.63302, 75.31743,
                 47.52159, 51.50512, 68.10810, -60.37985, 58.29240, 63.16882, 92.24195, 82.34067, 41.53038, 49.83299, 37.17257, 71.62839, 65.04276, 57.01429, 35.28390,
                 52.83621, 59.54155, 30.91311, 78.95270, 37.94634, 46.58579, 30.50557, 67.07761, 33.64414, -44.87616, 91.79820, 53.63460, 45.51329, 32.09087, 61.96482,
                 48.08096, 78.62163, 35.08889, 98.68567, 41.07923, 30.60649, 57.69076, 29.36967, 38.73110, 37.84127, 85.76051, 67.08190, 51.47189, 35.98709, -54.98000,
                 34.60094, 33.31931, 76.29741, 84.79908, 38.14804, 39.55889, 96.16675, 62.34194, 47.04409, -62.09536, 76.14043, 33.02666, 94.00564, 53.68852, 31.73132,
                 63.22703, 26.67201, 72.73931, -77.93963, 39.11328, 82.91384, 50.36465, 31.60422, -56.01568, 92.26648, 38.68745, 67.77602, 41.94134, 36.17816, 57.93474,
                 29.68316, 73.77379, 44.08081, 85.78152, 65.89571, 51.10218, 78.84303, -80.56622, -69.28504, -71.90990, 57.84604, 51.47220, 42.17414, -62.23792, 72.98064,
                 40.98654, 83.77524, 42.83898, 79.68969, -63.46598, -65.09480, 45.18274, 71.86100, -68.20659, -69.78288, -53.65702, 35.70765, -56.32207, -53.92707, 75.51740,
                 61.24951, -76.91300, -61.20788, 83.65608, -52.39311, 56.52574, 42.38686, 44.15033, -65.45875, 76.32325, 34.84570, -53.88628, 70.13359, 50.31694, -48.90561,
                 -71.56121, -62.40373, 44.82710, 61.05060, 40.77808, 81.06731, 52.35672, 33.41591, 34.32866, -64.06619, 39.38950, -48.12136, 51.74043, -52.23378, -69.81853,
                 78.23132, 53.79857, 72.25755, 40.73339, 60.65922, 43.76330, 47.55786, -46.33360, 38.63039, 66.36261, 39.85121, 41.35095, 54.76625, 80.35763, 61.44830,
                 43.58614, 45.35274, -43.61361, 50.88073, 72.04698, -49.60022, -56.37487, 48.96206, 42.82675, 48.90561, 59.53004, -39.93166, 85.53775, 41.55653, 56.46341,
                 79.78497, 50.17230, 67.60904, 44.77366, 47.45894, -62.74491, 57.73427, 61.68409, 74.56752, -44.63916, 58.97305, 49.70444, 51.26106, 57.01699, 86.40865,
                 73.61289, -65.20768, 63.42481, 55.80264, 74.97227, -80.45364, 64.28639, 83.00309, -67.30757, -50.91488, 60.46133, -76.38780, 70.86505, 97.19144, 81.65832,
                 67.42034, 66.74744, 77.86342, 88.41234, 68.08210, 84.76199, -57.80302, -59.69768, 98.95125, -71.39877, 53.28076, -90.53726, -72.39165, -68.50196, -72.27422,
                 -80.52008, -68.54955, -66.72530, -77.07195, -80.17271, -83.93850, -63.25588, -76.28748, -70.02824, 51.10953, -71.14010, 40.46037, -72.17862, -62.57266, 57.67305,
                 -55.90097, -55.39778, -66.67122, -56.42756, -61.20298, 39.14313, -74.94019, -65.25053, 50.15070, -67.58988, 63.21401, -57.96498, -62.43940, -75.79812, -51.70963,
                 -68.06333, -51.82955, -48.40215, 74.88232, 37.92051, -73.54655, 68.43147, 62.15394, -63.81999, -51.80399, -56.94498, -68.26764, -48.44406, -42.89581, 72.97069,
                 52.43054, 44.65180, -60.19297, -55.56035, -44.19917, 49.11250, 56.88895, -79.00633, -50.70163, -49.26376, 44.74917, 60.85039, 70.45090, -68.78874, -49.40788,
                 -59.14884, 32.55601, 54.20772, -47.36248, 33.61630, 50.39182, 41.55033, 62.29585, 67.79165, -68.78558, -53.99483, 42.54096, 46.58633, 70.58747, 58.49726,
                 -63.52971, 54.37175, 65.23877, -49.07492, 41.10413, 61.03113, 44.97764, -66.87476, -77.32278, 49.24589, 72.17169, -57.06608, -50.47158, -43.81847, -47.49533,
                 60.17817, 71.54758, 56.43004, -60.83717, 68.67465, 56.98453, -51.91766, 56.47203, 61.09518, -57.10227, -47.87032, 82.64731, 66.26489, 77.77807, 59.42314,
                 73.24015, -45.39070, -46.03222, 67.89339, 80.60991, 87.66147, 66.61122, 77.62832, 64.50831, 76.77484, -50.84583, 72.57423, -54.62405, 74.23006, 78.38152,
                 73.46824, 86.50457), 
               tolerance = 0.00001)

  # works on cell types when present and no negative distances returned
  expect_message(cur_spe_2 <- distToCells(object = cur_spe,
                                          x_cells = cur_spe$CellType == "celltype_B",
                                          coords = c("Pos_X","Pos_Y"),
                                          metric = "median",
                                          img_id = "ImageName",
                                          return_neg = FALSE), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_true(is(cur_spe_2, "SingleCellExperiment"))
  expect_s4_class(cur_spe_2 , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_spe_2)))
  expect_true(class(cur_spe_2$distToCells) == "numeric")
  expect_true(min(cur_spe_2$distToCells) > 0)

  expect_equal(cur_spe_1[,!cur_spe$CellType == "celltype_B"]$distToCells,cur_spe_2[,!cur_spe$CellType == "celltype_B"]$distToCells)

  expect_equal(sum(cur_spe_1$distToCells < 0),sum(cur_spe$CellType == "celltype_B"))

  # compare results from SingleCellExperiment and SpatialExperiment
  expect_equal(cur_sce$distToCells,cur_spe_1$distToCells)

  expect_equal(cur_sce_2$distToCells,cur_spe_2$distToCells)

  # Works when all cells of an image belong to one batch
  expect_message(cur_sce <- distToCells(object = pancreasSCE,
                                        x_cells = pancreasSCE$ImageName == "J02_imc.tiff",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "median",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")

  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(all(is.na(cur_sce$distToCells)))

  cur_sce$CellType[cur_sce$ImageName == "J02_imc.tiff"] <- "celltype_A"
  expect_message(cur_sce <- distToCells(object = cur_sce,
                                        x_cells = cur_sce$CellType == "celltype_A",
                                        coords = c("Pos_X","Pos_Y"),
                                        metric = "median",
                                        img_id = "ImageName"), regexp = "The returned object is ordered by the 'ImageName' entry.")
  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_true("distToCells" %in% names(colData(cur_sce)))
  expect_true(all(is.na(cur_sce$distToCells[cur_sce$ImageName == "J02_imc.tiff"])))
  expect_true(all(!is.na(cur_sce$distToCells[cur_sce$ImageName != "J02_imc.tiff"])))


  # Error
  expect_error(cur_sce_4 <- distToCells(object = pancreasSCE[,pancreasSCE$ImageName == "J02_imc.tiff"],
                                           x_cells = pancreasSCE$CellType == "celltype_A",
                                           coords = c("Pos_X","Pos_Y"),
                                           img_id = "ImageName",
                                           return_neg = FALSE),
               regexp = "Length of 'x_cells' must match the number of cells in 'object'.")

  expect_error(distToCells(object = "test"),
               regexp = "'object' not of type 'SingleCellExperiment'.",
               fixed = TRUE)
  expect_error(distToCells(object = pancreasSCE[,pancreasSCE$ImageName == "test"], x_cells = pancreasSCE[,pancreasSCE$ImageName == "test"]$CellType ==  "celltype_B",name = "test",coords = c("Pos_X","Pos_Y"),
               img_id = "ImageName",return_neg = TRUE),
               regexp = "'object' must contain at least one cell",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = "test"),
               regexp = "'x_cells' must all be logical.",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = TRUE),
               regexp = "Length of 'x_cells' must match the number of cells in 'object'.",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = pancreasSCE$CellType ==  "celltype_B",name = TRUE),
               regexp = "'name' must be a single string.",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = pancreasSCE$CellType ==  "celltype_B",name = 1),
               regexp = "'name' must be a single string.",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = pancreasSCE$CellType ==  "celltype_B",metric = "test"),
               regexp = "'metric' not supported. Must be one of 'min', 'max', 'mean' or 'median'",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = pancreasSCE$CellType ==  "celltype_B",name = "test",coords = c(1,2)),
               regexp = "'coords' must be a character vector of length 2.",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = pancreasSCE$CellType ==  "celltype_B",name = "test",coords = c("A","B")),
               regexp = "'coords' not in colData(object).",
               fixed = TRUE)
  expect_error(distToCells(cur_spe, x_cells = pancreasSCE$CellType ==  "celltype_B",name = "test",coords = c("A","B")),
               regexp = "'coords' not in spatialCoords(object).",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = pancreasSCE$CellType ==  "celltype_B",name = "test",coords = c("Pos_X","Pos_Y"),img_id = 1),
               regexp = "'img_id' must be a single string.",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = pancreasSCE$CellType ==  "celltype_B",name = "test",coords = c("Pos_X","Pos_Y"),img_id = "test"),
               regexp = "'img_id' not in colData(object).",
               fixed = TRUE)
  expect_error(distToCells(pancreasSCE, x_cells = pancreasSCE$CellType ==  "celltype_B",name = "test",coords = c("Pos_X","Pos_Y"),
                              img_id = "ImageName",return_neg = 1),
               regexp = "'return_neg' is not of type logical.",
               fixed = TRUE)
  expect_error(distToCells(cur_spe, x_cells = pancreasSCE$CellType ==  "celltype_B",name = "test",coords = c("Pos_1","Pos_Y")),
               regexp = "'coords' not in spatialCoords(object).",
               fixed = TRUE)
})


