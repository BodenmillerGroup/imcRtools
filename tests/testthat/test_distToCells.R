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
  
  # Check numerical values
  expect_equal(cur_sce$distToCells %>% unname, c(79.788514, 20.857982, 53.368632, 62.713162, 19.560774, 41.448102, 62.462909, 68.450475, 31.699103, 54.291193, 10.058905, 46.570420, 10.015750,
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
                                      15.379208, 24.752166, 9.317498, 19.900208, -11.012489, 23.229380, -9.911600, 11.012489, 29.874798, 13.596656, 29.695571)
               , 
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
  # TODO: for all cases, enforces that a cells not of interest and a cell of interect do not have the same exact coordinates

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
  # TODO: for all cases, enforces that a cells not of interest and a cell of interect do not have the same exact coordinates

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
  # TODO: for all cases, enforces that a cells not of interest and a cell of interect do not have the same exact coordinates

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


