test_that("readSCEfromTIFF function reads in correct objects.", {
  path <- system.file("extdata/spillover_tiff/img", package = "imcRtools")
  image_df_path <- system.file("extdata/spillover_tiff/images.csv", package = "imcRtools")
  panel_df_path <- system.file("extdata/spillover_tiff/panel.csv", package = "imcRtools")
  
  # Read in .tiff
  expect_silent(cur_sce <- readSCEfromTIFF(path, image_df_path, panel_df_path, verbose = FALSE))
  expect_equal(rowData(cur_sce)$channel_name, c("Dy161Di", "Dy162Di", 
                                                "Dy163Di","Dy164Di"))
  expect_equal(rowData(cur_sce)$marker_name, c("Dy161", "Dy162", 
                                               "Dy163","Dy164"))
  expect_equal(cur_sce$sample_id, rep(c("Dy161", "Dy162", 
                                        "Dy163","Dy164"), each = 500))
  expect_equal(cur_sce$sample_metal, rep(c("Dy", "Dy", 
                                           "Dy","Dy"), each = 500))
  expect_equal(cur_sce$sample_mass, rep(c("161", "162", 
                                          "163","164"), each = 500))
  expect_equal(assayNames(cur_sce), "counts")
  expect_equal(dim(cur_sce), c(4, 2000))
  
  # Verbose output
  cur_out <- capture_output(cur_sce <- readSCEfromTIFF(path, image_df_path, panel_df_path))
  expect_equal(cur_out, "Spotted channels:  Dy161, Dy162, Dy163, Dy164\nAcquired channels:  Dy161, Dy162, Dy163, Dy164\nChannels spotted but not acquired:  \nChannels acquired but not spotted:  ")
  
  # Error
  expect_error(readSCEfromTIFF("test", image_df_path, panel_df_path, verbose = FALSE), 
               regexp = "Image folder path does not exist.",
               fixed = TRUE)
  
  expect_error(readSCEfromTIFF(path, "test", panel_df_path, verbose = FALSE), 
               regexp = "Image dataframe path does not exist.",
               fixed = TRUE)
  
  expect_error(readSCEfromTIFF(path, image_df_path, "test", verbose = FALSE), 
               regexp = "Panel dataframe path does not exist.",
               fixed = TRUE)
})
