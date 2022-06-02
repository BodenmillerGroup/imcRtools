test_that("read_steinbock function works", {
    path <- system.file("extdata/mockData/steinbock", package = "imcRtools")
  
    # SpatialExperiment
    cur_spe <- read_steinbock(path)
    
    expect_equal(length(int_metadata(cur_spe)), 3)
    
    expect_s4_class(cur_spe, "SpatialExperiment")
    
    expect_equal(rownames(cur_spe), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 218))
    expect_equal(names(rowData(cur_spe)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "area", 
                                            "major_axis_length", 
                                            "minor_axis_length", "eccentricity",
                                            "width_px", "height_px"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv, show_col_types = FALSE)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_spe), t(cur_counts[,-1]))
    expect_equal(counts(cur_spe)[1:10], c(0.065453581, 0.046153846, 0.076923077, 
                                          0.096965075, 0.750148960, 0.084382993, 
                                          0.058186264, 0.004239203, 0.055236769, 
                                          0.658203733), 
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv, show_col_types = FALSE)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_spe$area, cur_morph$area)
    expect_equal(cur_spe$major_axis_length, cur_morph$major_axis_length)
    expect_equal(cur_spe$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_spe$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), cur_morph$`centroid-1`)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), cur_morph$`centroid-0`)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"), show_col_types = FALSE)
    expect_equal(rowData(cur_spe)$name, cur_panel$name)
    expect_equal(rowData(cur_spe)$channel, cur_panel$channel)
    expect_equal(rowData(cur_spe)$keep, cur_panel$keep)
    expect_equal(rowData(cur_spe)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(colPairNames(cur_spe), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData1_3.csv"), show_col_types = FALSE)
    
    expect_equal(cur_test$Object + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData1_3")])
    expect_equal(cur_test$Neighbor + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData1_3")])
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData2_1.csv"), show_col_types = FALSE)
    
    expect_equal(cur_test$Object + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData2_1")])
    expect_equal(cur_test$Neighbor + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData2_1")])
    
    # SingleCellExperiment
    cur_sce <- read_steinbock(path, return_as = "sce")
    
    expect_s4_class(cur_sce, "SingleCellExperiment")
    
    expect_equal(rownames(cur_sce), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 218))
    expect_equal(names(rowData(cur_sce)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", 
                                            "Pos_Y", "area", 
                                            "major_axis_length", 
                                            "minor_axis_length", "eccentricity",
                                            "width_px", "height_px"))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv, show_col_types = FALSE)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_sce), t(cur_counts[,-1]))
    expect_equal(counts(cur_sce)[1:10], c(0.065453581, 0.046153846, 0.076923077, 
                                          0.096965075, 0.750148960, 0.084382993, 
                                          0.058186264, 0.004239203, 0.055236769, 
                                          0.658203733), 
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv, show_col_types = FALSE)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_sce$area, cur_morph$area)
    expect_equal(cur_sce$major_axis_length, cur_morph$major_axis_length)
    expect_equal(cur_sce$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_sce$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(cur_sce$Pos_X), cur_morph$`centroid-1`)
    expect_equal(as.numeric(cur_sce$Pos_Y), cur_morph$`centroid-0`)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"), show_col_types = FALSE)
    expect_equal(rowData(cur_sce)$name, cur_panel$name)
    expect_equal(rowData(cur_sce)$channel, cur_panel$channel)
    expect_equal(rowData(cur_sce)$keep, cur_panel$keep)
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_sce)$deepcell, cur_panel$deepcell)    
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(colPairNames(cur_sce), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData1_3.csv"), show_col_types = FALSE)
    
    expect_equal(cur_test$Object + sum(cur_sce$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData1_3")])
    expect_equal(cur_test$Neighbor + sum(cur_sce$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData1_3")])
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData2_1.csv"), show_col_types = FALSE)
    
    expect_equal(cur_test$Object + sum(cur_sce$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData2_1")])
    expect_equal(cur_test$Neighbor + sum(cur_sce$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData2_1")])
    
    # Test other import settings
    cur_spe <- read_steinbock(path, regionprops_folder = NULL)
    
    expect_equal(rownames(cur_spe), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 218))
    expect_equal(names(rowData(cur_spe)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", 
                                            "width_px", "height_px"))
    expect_null(spatialCoordsNames(cur_spe))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv, show_col_types = FALSE)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_spe), t(cur_counts[,-1]))
    expect_equal(counts(cur_spe)[1:10], c(0.065453581, 0.046153846, 0.076923077, 
                                          0.096965075, 0.750148960, 0.084382993, 
                                          0.058186264, 0.004239203, 0.055236769, 
                                          0.658203733), 
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv, show_col_types = FALSE)
    cur_morph <- do.call("rbind", cur_morph)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"), show_col_types = FALSE)
    expect_equal(rowData(cur_spe)$name, cur_panel$name)
    expect_equal(rowData(cur_spe)$channel, cur_panel$channel)
    expect_equal(rowData(cur_spe)$keep, cur_panel$keep)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(colPairNames(cur_sce), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData1_3.csv"), show_col_types = FALSE)
    
    expect_equal(cur_test$Object + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData1_3")])
    expect_equal(cur_test$Neighbor + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData1_3")])
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData2_1.csv"), show_col_types = FALSE)
    
    expect_equal(cur_test$Object + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData2_1")])
    expect_equal(cur_test$Neighbor + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData2_1")])
    
    cur_sce <- read_steinbock(path, return_as = "sce", regionprops_folder = NULL)
    
    expect_equal(rownames(cur_sce), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 218))
    expect_equal(names(rowData(cur_sce)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber",
                                            "width_px", "height_px"))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv, show_col_types = FALSE)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_sce), t(cur_counts[,-1]))
    expect_equal(counts(cur_sce)[1:10], c(0.065453581, 0.046153846, 0.076923077, 
                                          0.096965075, 0.750148960, 0.084382993, 
                                          0.058186264, 0.004239203, 0.055236769, 
                                          0.658203733), 
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv, show_col_types = FALSE)
    cur_morph <- do.call("rbind", cur_morph)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"), show_col_types = FALSE)
    expect_equal(rowData(cur_sce)$name, cur_panel$name)
    expect_equal(rowData(cur_sce)$channel, cur_panel$channel)
    expect_equal(rowData(cur_sce)$keep, cur_panel$keep)
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_sce)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(colPairNames(cur_sce), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData1_3.csv"), show_col_types = FALSE)
    
    expect_equal(cur_test$Object + sum(cur_sce$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData1_3")])
    expect_equal(cur_test$Neighbor + sum(cur_sce$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData1_3")])
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData2_1.csv"))
    
    expect_equal(cur_test$Object + sum(cur_sce$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData2_1")])
    expect_equal(cur_test$Neighbor + sum(cur_sce$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData2_1")])
    
    
    cur_spe <- read_steinbock(path, graphs_folder = NULL)
    
    expect_equal(rownames(cur_spe), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 218))
    expect_equal(names(rowData(cur_spe)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "area", 
                                            "major_axis_length", 
                                            "minor_axis_length", "eccentricity",
                                            "width_px", "height_px"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_spe), t(cur_counts[,-1]))
    expect_equal(counts(cur_spe)[1:10], c(0.065453581, 0.046153846, 0.076923077, 
                                          0.096965075, 0.750148960, 0.084382993, 
                                          0.058186264, 0.004239203, 0.055236769, 
                                          0.658203733), 
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_spe$area, cur_morph$area)
    expect_equal(cur_spe$major_axis_length, cur_morph$major_axis_length)
    expect_equal(cur_spe$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_spe$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), cur_morph$`centroid-1`)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), cur_morph$`centroid-0`)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"))
    expect_equal(rowData(cur_spe)$name, cur_panel$name)
    expect_equal(rowData(cur_spe)$channel, cur_panel$channel)
    expect_equal(rowData(cur_spe)$keep, cur_panel$keep)
    expect_equal(rowData(cur_spe)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(length(colPairs(cur_spe)), 0)
    
    cur_sce <- read_steinbock(path, return_as = "sce", graphs_folder = NULL)
    
    expect_s4_class(cur_sce, "SingleCellExperiment")
    
    expect_equal(rownames(cur_sce), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 218))
    expect_equal(names(rowData(cur_sce)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", 
                                            "Pos_Y", "area", 
                                            "major_axis_length", 
                                            "minor_axis_length", "eccentricity",
                                            "width_px", "height_px"))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_sce), t(cur_counts[,-1]))
    expect_equal(counts(cur_sce)[1:10], c(0.065453581, 0.046153846, 0.076923077, 
                                          0.096965075, 0.750148960, 0.084382993, 
                                          0.058186264, 0.004239203, 0.055236769, 
                                          0.658203733),
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_sce$area, cur_morph$area)
    expect_equal(cur_sce$major_axis_length, cur_morph$major_axis_length)
    expect_equal(cur_sce$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_sce$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(cur_sce$Pos_X), cur_morph$`centroid-1`)
    expect_equal(as.numeric(cur_sce$Pos_Y), cur_morph$`centroid-0`)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"))
    expect_equal(rowData(cur_sce)$name, cur_panel$name)
    expect_equal(rowData(cur_sce)$channel, cur_panel$channel)
    expect_equal(rowData(cur_sce)$keep, cur_panel$keep)
    expect_equal(rowData(cur_sce)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(length(colPairs(cur_sce)), 0)
    
    cur_spe <- read_steinbock(path, graphs_folder = NULL, regionprops_folder = NULL)
    
    expect_equal(rownames(cur_spe), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 218))
    expect_equal(names(rowData(cur_spe)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber",
                                            "width_px", "height_px"))
    expect_null(spatialCoordsNames(cur_spe))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_spe), t(cur_counts[,-1]))
    expect_equal(counts(cur_spe)[1:10], c(0.065453581, 0.046153846, 0.076923077, 
                                          0.096965075, 0.750148960, 0.084382993, 
                                          0.058186264, 0.004239203, 0.055236769, 
                                          0.658203733),
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"))
    expect_equal(rowData(cur_spe)$name, cur_panel$name)
    expect_equal(rowData(cur_spe)$channel, cur_panel$channel)
    expect_equal(rowData(cur_spe)$keep, cur_panel$keep)
    expect_equal(rowData(cur_spe)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(length(colPairs(cur_spe)), 0)
    
    cur_sce <- read_steinbock(path, return_as = "sce", graphs_folder = NULL, regionprops_folder = NULL)
    
    expect_s4_class(cur_sce, "SingleCellExperiment")
    
    expect_equal(rownames(cur_sce), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 218))
    expect_equal(names(rowData(cur_sce)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber",
                                            "width_px", "height_px"))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_sce), t(cur_counts[,-1]))
    expect_equal(counts(cur_sce)[1:10], c(0.065453581, 0.046153846, 0.076923077, 
                                          0.096965075, 0.750148960, 0.084382993, 
                                          0.058186264, 0.004239203, 0.055236769, 
                                          0.658203733),
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"))
    expect_equal(rowData(cur_sce)$name, cur_panel$name)
    expect_equal(rowData(cur_sce)$channel, cur_panel$channel)
    expect_equal(rowData(cur_sce)$keep, cur_panel$keep)
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_sce)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(length(colPairs(cur_sce)), 0)
    
    cur_spe <- read_steinbock(path, pattern = "mockData2")
    
    expect_s4_class(cur_spe, "SpatialExperiment")
    
    expect_equal(rownames(cur_spe), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 35))
    expect_equal(names(rowData(cur_spe)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "area", 
                                            "major_axis_length", 
                                            "minor_axis_length", "eccentricity",
                                            "width_px", "height_px"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE, pattern = "mockData2")
    cur_counts <- lapply(cur_files, readr::read_csv)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_spe), t(cur_counts[,-1]))
    expect_equal(counts(cur_spe)[1:10], c(0.38237867, 0.04112295, 6.88474564, 
                                          5.65280831, 2.44857732, 0.53002120, 
                                          0.08389935, 2.76029619, 8.47630586, 
                                          3.58901901),
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE, pattern = "mockData2")
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_spe$area, cur_morph$area)
    expect_equal(cur_spe$major_axis_length, cur_morph$major_axis_length)
    expect_equal(cur_spe$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_spe$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), cur_morph$`centroid-1`)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), cur_morph$`centroid-0`)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"))
    expect_equal(rowData(cur_spe)$name, cur_panel$name)
    expect_equal(rowData(cur_spe)$channel, cur_panel$channel)
    expect_equal(rowData(cur_spe)$keep, cur_panel$keep)
    expect_equal(rowData(cur_spe)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(colPairNames(cur_spe), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData2_3.csv"))
    
    expect_equal(cur_test$Object + sum(cur_spe$sample_id %in% c("20210305_NE_mockData2_1", "20210305_NE_mockData2_2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData2_3")])
    expect_equal(cur_test$Neighbor + sum(cur_spe$sample_id %in% c("20210305_NE_mockData2_1", "20210305_NE_mockData2_2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData2_3")])
    
    cur_sce <- read_steinbock(path, pattern = "mockData2", return_as = "sce")
    
    expect_s4_class(cur_sce, "SingleCellExperiment")
    
    expect_equal(rownames(cur_sce), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 35))
    expect_equal(names(rowData(cur_sce)), 
                 c("channel", "name", "keep", "ilastik", "deepcell", "Tube.Number"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", 
                                            "Pos_Y", "area", 
                                            "major_axis_length", 
                                            "minor_axis_length", "eccentricity",
                                            "width_px", "height_px"))
    
    cur_files <- list.files(file.path(path, "intensities"), full.names = TRUE, pattern = "mockData2")
    cur_counts <- lapply(cur_files, readr::read_csv)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_sce), t(cur_counts[,-1]))
    expect_equal(counts(cur_sce)[1:10], c(0.38237867, 0.04112295, 6.88474564, 
                                          5.65280831, 2.44857732, 0.53002120, 
                                          0.08389935, 2.76029619, 8.47630586, 
                                          3.58901901),
                 tolerance = 10e-6)
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE, pattern = "mockData2")
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_sce$area, cur_morph$area)
    expect_equal(cur_sce$major_axis_length, cur_morph$major_axis_length)
    expect_equal(cur_sce$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_sce$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(cur_sce$Pos_X), cur_morph$`centroid-1`)
    expect_equal(as.numeric(cur_sce$Pos_Y), cur_morph$`centroid-0`)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"))
    expect_equal(rowData(cur_sce)$name, cur_panel$name)
    expect_equal(rowData(cur_sce)$channel, cur_panel$channel)
    expect_equal(rowData(cur_sce)$keep, cur_panel$keep)
    expect_equal(rowData(cur_sce)$deepcell, cur_panel$deepcell)
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_equal(colPairNames(cur_sce), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- readr::read_csv(file.path(path, "neighbors", "20210305_NE_mockData2_3.csv"))
    
    expect_equal(cur_test$Object + sum(cur_sce$sample_id %in% c("20210305_NE_mockData2_1", "20210305_NE_mockData2_2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData2_3")])
    expect_equal(cur_test$Neighbor + sum(cur_sce$sample_id %in% c("20210305_NE_mockData2_1", "20210305_NE_mockData2_2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "20210305_NE_mockData2_3")])
    
    cur_spe <- read_steinbock(path, panel = NULL)
    
    expect_equal(length(rowData(cur_spe)), 0)
    
    cur_sce <- read_steinbock(path, panel = NULL, return_as = "sce")
    
    expect_equal(length(rowData(cur_sce)), 0)
    
    cur_spe <- read_steinbock(path, extract_coords_from = c("area", "major_axis_length"))
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "centroid.0",
                                            "centroid.1", "minor_axis_length", "eccentricity",
                                            "width_px", "height_px"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_spe$`centroid.1`, cur_morph$`centroid-1`)
    expect_equal(cur_spe$`centroid.0`, cur_morph$`centroid-0`)
    expect_equal(cur_spe$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_spe$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), cur_morph$area)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), cur_morph$major_axis_length)
    
    cur_sce <- read_steinbock(path, return_as = "sce", extract_coords_from = c("area", "major_axis_length"))
    
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", "Pos_Y", "centroid.0",
                                            "centroid.1", "minor_axis_length", "eccentricity",
                                            "width_px", "height_px"))
    
    cur_files <- list.files(file.path(path, "regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_sce$`centroid.0`, cur_morph$`centroid-0`)
    expect_equal(cur_sce$`centroid.1`, cur_morph$`centroid-1`)
    expect_equal(cur_sce$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_sce$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(cur_sce$Pos_X), cur_morph$area)
    expect_equal(as.numeric(cur_sce$Pos_Y), cur_morph$major_axis_length)
    
    cur_spe <- read_steinbock(path, extract_names_from = "channel")
    
    expect_true(all(is.na(rowData(cur_spe)["Laminin",])))
    expect_equal(as.character(as.matrix(rowData(cur_spe)["Ag107",])), 
                 c("Ag107", "Ag107", "1", "1", NA, NA))  
    
    cur_spe <- read_steinbock(path, image_file = NULL)
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "area", 
                                            "major_axis_length", "minor_axis_length", "eccentricity"))
    
    cur_sce <- read_steinbock(path, image_file = NULL, return_as = "sce")
    
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", "Pos_Y", "area", 
                                            "major_axis_length", "minor_axis_length", "eccentricity"))
    
    cur_spe <- read_steinbock(path)
    cur_images_file <- readr::read_csv(file.path(path, "images.csv"), show_col_types = FALSE)
    
    cur_df <- unique(colData(cur_spe)[,c("sample_id", "width_px", "height_px")])
    expect_equal(cur_images_file$image, paste0(cur_df$sample_id, ".tiff"))
    expect_equal(cur_images_file$width_px, cur_df$width_px)
    expect_equal(cur_images_file$height_px, cur_df$height_px)
    
    cur_spe <- read_steinbock(path, extract_imagemetadata_from = c("recovered", "acquisition_description"))
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "area", 
                                            "major_axis_length", "minor_axis_length", 
                                            "eccentricity", "recovered", "acquisition_description"))
    
    cur_df <- unique(colData(cur_spe)[,c("sample_id", "recovered", "acquisition_description")])
    expect_equal(cur_images_file$image, paste0(cur_df$sample_id, ".tiff"))
    expect_equal(cur_images_file$recovered, cur_df$recovered)
    expect_equal(cur_images_file$acquisition_description, cur_df$acquisition_description)
    
    
    
    # Parallelisation
    #cur_spe <- read_steinbock(path, BPPARAM = BiocParallel::bpparam())
    
    # Error
    expect_error(cur_spe <- read_steinbock(path = "test"),
                 "'path' doesn't exist.", 
                 fixed = TRUE) 
    expect_error(cur_spe <- read_steinbock(path = c("test", "test2")),
                 "'path' must be a single string.", 
                 fixed = TRUE) 
    expect_error(cur_spe <- read_steinbock(path = 1),
                 "'path' must be a single string.", 
                 fixed = TRUE) 
    
    expect_error(cur_spe <- read_steinbock(path, intensities_folder = NULL),
                 "'intensities_folder' must be specified.", 
                 fixed = TRUE)
    expect_error(cur_spe <- read_steinbock(path, intensities_folder = "test"),
                 "'intensities_folder' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, intensities_folder = c("test", "test2")),
                 "'intensities_folder' must be a single string.", 
                 fixed = TRUE)
    expect_error(cur_spe <- read_steinbock(path, intensities_folder = 1),
                 "'intensities_folder' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, regionprops_folder = "test"),
                 "'regionprops_folder' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, regionprops_folder = c("test", "test2")),
                 "'regionprops_folder' must be a single string.", 
                 fixed = TRUE)
    expect_error(cur_spe <- read_steinbock(path, regionprops_folder = 1),
                 "'regionprops_folder' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, graphs_folder = "test"),
                 "'graphs_folder' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, graphs_folder = c("test", "test2")),
                 "'graphs_folder' must be a single string.", 
                 fixed = TRUE)
    expect_error(cur_spe <- read_steinbock(path, graphs_folder = 1),
                 "'graphs_folder' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, pattern = "test"),
                 "No files were read in.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, extract_cellid_from = "test"),
                 "'extract_cellid_from' not in intensities files.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, extract_cellid_from = c("test", "test2")),
                 "'extract_cellid_from' must be a single string.", 
                 fixed = TRUE)
    expect_error(cur_spe <- read_steinbock(path, extract_cellid_from = 1),
                 "'extract_cellid_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, extract_coords_from =  "test"),
                 "'coords' not in regionprops files.", 
                 fixed = TRUE)
    expect_silent(cur_spe <- read_steinbock(path, extract_coords_from =  "test", regionprops_folder = NULL))
    
    expect_error(cur_spe <- read_steinbock(path, extract_coords_from =  1),
                 "'extract_coords_from' must be characters.", 
                 fixed = TRUE)
    
    expect_warning(cur_spe <- read_steinbock(path, panel = "test"),
                   "'panel_file' does not exist.", 
                   fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, panel = c("test", "test2")),
                 "'panel_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, panel = 1),
                 "'panel_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, extract_names_from = "test"),
                 "'extract_names_from' not in panel file.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, extract_names_from = c("test", "test2")),
                 "'extract_names_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, extract_names_from = 1),
                 "'extract_names_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, image_file = 1),
                 "'image_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, image_file = "test"),
                 "'image_file' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, extract_imagemetadata_from = 1),
                 "'extract_imagemetadata_from' should only contain characters.", 
                 fixed = TRUE)
    
    expect_error(cur_spe <- read_steinbock(path, extract_imagemetadata_from = c(1, "test")),
                 "'extract_imagemetadata_from' not in images file.", 
                 fixed = TRUE)
})

test_that("read_steinbock function works when files are missing", {
    skip_on_os(os = "windows")
    path <- system.file("extdata/mockData/steinbock", package = "imcRtools")
    
    # Move files to tmp location
    cur_path <- tempdir()
    file.copy(path, cur_path, recursive = TRUE)
    
    # Remove regionprobs folder
    file.remove(list.files(paste0(cur_path, "/steinbock/regionprops"), 
                           full.names = TRUE))
    
    cur_spe <- read_steinbock(paste0(cur_path, "/steinbock/"))
    
    expect_s4_class(cur_spe, "SpatialExperiment")
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber",
                                            "width_px", "height_px"))
    
    expect_equal(colPairNames(cur_spe), "neighborhood")

    # Remove graphs folder
    file.remove(list.files(paste0(cur_path, "/steinbock/neighbors"), 
                           full.names = TRUE))
    
    cur_spe <- read_steinbock(paste0(cur_path, "/steinbock/"))
    
    expect_s4_class(cur_spe, "SpatialExperiment")
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber",
                                            "width_px", "height_px"))
    
    expect_error(colPair(cur_spe), 
                 regex = "no available entries for 'colPair(<SpatialExperiment>, ...)'",
                 fixed = TRUE)
    
    # Copy panel
    file.copy(paste0(cur_path, "/steinbock/panel.csv"), 
              paste0(cur_path, "/steinbock/panel_2.csv"))
    
    expect_silent(cur_spe <- read_steinbock(paste0(cur_path, "/steinbock/"), 
                              panel_file = paste0(cur_path, "/steinbock/panel_2.csv")))

})
