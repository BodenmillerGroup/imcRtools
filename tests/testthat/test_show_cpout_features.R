test_that("show_cpout_features function works", {
    path <- system.file("extdata/mockData/cpout", package = "imcRtools")
    
    # Works
    expect_silent(show_cpout_features(path))
    expect_silent(cur_tab <- show_cpout_features(path))
    expect_s3_class(cur_tab, "datatables")
    
    cur_df <- cur_tab$x$data
    expect_equal(colnames(cur_df), c(" ", "column_name", "category", "image_name", "object_name", "feature_name", "channel", 
                                     "parameters", "channel_id", "data_type"))
    expect_equal(unique(cur_df$channel_id), c(NA, "Ag107", "Pr141", "Sm147", "Eu153", "Yb172", "CellCenter", 
                                              "CellBorder", "Background"))
    expect_equal(head(cur_df$column_name), c("Location_Center_X", "Location_Center_Y", "Location_Center_Z",
                                             "Number_Object_Number", "Neighbors_NumberOfNeighbors_8",
                                             "Neighbors_PercentTouching_8"))
    
    expect_silent(show_cpout_features(path, display = "image_features"))
    expect_silent(cur_tab <- show_cpout_features(path, display = "image_features"))
    expect_s3_class(cur_tab, "datatables")
    
    cur_df <- cur_tab$x$data
    expect_equal(colnames(cur_df), c(" ", "column_name", "category", "image_name", "object_name", "feature_name", "channel", 
                                     "parameters", "channel_id", "data_type"))
    expect_equal(unique(cur_df$channel_id), c(NA, "Ag107", "Pr141", "Sm147", "Eu153", "Yb172", "CellCenter", 
                                              "CellBorder", "Background"))
    expect_equal(head(cur_df$column_name), c("Group_Number", "Group_Index", "ModuleError_01Images", "ExecutionTime_01Images", 
                                             "Metadata_", "Metadata_template"))
    
    
    # Error
    expect_error(show_cpout_features("test"),
                 regexp = "'var_cell.csv' does not exist in test",
                 fixed = TRUE)
    expect_error(show_cpout_features(path, cell_features = "test"),
                 regexp = paste0("'test' does not exist in ", path),
                 fixed = TRUE)
    expect_error(show_cpout_features("test", display = "image_features"),
                 regexp = "'var_Image.csv' does not exist in test",
                 fixed = TRUE)
    expect_error(show_cpout_features(path, image_features = "test", display = "image_features"),
                 regexp = paste0("'test' does not exist in ", path),
                 fixed = TRUE)
    expect_error(show_cpout_features(path, display = "test"))
})
