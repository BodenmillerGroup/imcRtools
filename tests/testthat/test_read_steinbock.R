test_that("read_steinbock function works", {
    path <- system.file("extdata/mockData/steinbock", package = "imcRtools")
  
    # SpatialExperiment
    cur_spe <- read_steinbock(path)
    
    expect_equal(rownames(cur_spe), c("Ag107", "Cytokeratin 5", "Laminin", 
                                      "YBX1", "H3K27Ac"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 348))
    expect_equal(names(rowData(cur_spe)), 
                 c("channel", "name", "keep", "ilastik", "Tube.Number"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectId", "area", 
                                            "major_axis_length", 
                                            "minor_axis_length", "eccentricity"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    cur_files <- list.files(file.path(path, "cell_intensities"), full.names = TRUE)
    cur_counts <- lapply(cur_files, readr::read_csv)
    cur_counts <- do.call("rbind", cur_counts)
    
    expect_equal(counts(cur_spe), t(cur_counts[,-1]))
    expect_equal(counts(cur_spe)[1:10], c(0.19277269, 0.10526316, 0.10526316, 0.33692434, 0.53972906,
                                          0.05992229, 0.05633803, 0.04225352, 0.09859155, 0.72663587))
    
    cur_files <- list.files(file.path(path, "cell_regionprops"), full.names = TRUE)
    cur_morph <- lapply(cur_files, readr::read_csv)
    cur_morph <- do.call("rbind", cur_morph)
    
    expect_equal(cur_spe$area, cur_morph$area)
    expect_equal(cur_spe$major_axis_length, cur_morph$major_axis_length)
    expect_equal(cur_spe$minor_axis_length, cur_morph$minor_axis_length)
    expect_equal(cur_spe$eccentricity, cur_morph$eccentricity)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), cur_morph$`centroid-0`)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), cur_morph$`centroid-1`)
    
    cur_panel <- readr::read_csv(file.path(path, "panel.csv"))
    expect_equal(rowData(cur_spe)$name, cur_panel$name)
    expect_equal(rowData(cur_spe)$channel, cur_panel$channel)
    expect_equal(rowData(cur_spe)$keep, cur_panel$keep)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighbourhood"))
    cur_test <- readr::read_csv(file.path(path, "cell_graphs", "20210305_NE_mockData1_3.csv"))
    
    expect_equal(cur_test$Object1 + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData1_3")])
    expect_equal(cur_test$Object2 + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData1_3")])
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighbourhood"))
    cur_test <- readr::read_csv(file.path(path, "cell_graphs", "20210305_NE_mockData2_1.csv"))
    
    expect_equal(cur_test$Object1 + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData2_1")])
    expect_equal(cur_test$Object2 + sum(cur_spe$sample_id %in% c("20210305_NE_mockData1_1", "20210305_NE_mockData1_2", "20210305_NE_mockData1_3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "20210305_NE_mockData2_1")])
    
    # SingleCellExperiment
    
})
