test_that("read_steinbock function works", {
    path <- system.file("extdata/mockData/steinbock", package = "imcRtools")
  
    # SpatialExperiment
    cur_spe <- read_steinbock(path)
    
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
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), cur_morph$`centroid-0`)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), cur_morph$`centroid-1`)
    
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
})
