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


