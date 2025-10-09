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
                                        "Dy163","Dy164"), each = 20000))
  expect_equal(cur_sce$sample_metal, rep(c("Dy", "Dy", 
                                           "Dy","Dy"), each = 20000))
  expect_equal(cur_sce$sample_mass, rep(c("161", "162", 
                                          "163","164"), each = 20000))
  expect_equal(assayNames(cur_sce), "counts")
  expect_equal(dim(cur_sce), c(4, 80000))
  
  # Verbose output
  cur_out <- capture_output(cur_sce <- readSCEfromTIFF(path, image_df_path, panel_df_path))
  expect_equal(cur_out, "Spotted channels:  Dy161, Dy162, Dy163, Dy164\nAcquired channels:  Dy161, Dy162, Dy163, Dy164\nChannels spotted but not acquired:  \nChannels acquired but not spotted:  ")
  
  # Compare spillover matrix created with readSCEfromTXT to readSCEfromTIFF
  path_txt <- system.file("extdata/spillover_tiff/txt", package = "imcRtools")
  txt_sce <- readSCEfromTXT(path_txt, verbose = FALSE)
  txt_sce <- readSCEfromTXT(path_txt, verbose = FALSE)
  assay(txt_sce, "exprs") <- asinh(counts(txt_sce) / 5)
  bc_key_txt <- as.numeric(unique(txt_sce$sample_mass))
  bc_key_txt <- bc_key_txt[order(bc_key_txt)]
  txt_sce <- assignPrelim(txt_sce, bc_key = bc_key_txt)
  txt_sce <- estCutoffs(txt_sce)
  txt_sce <- applyCutoffs(txt_sce)
  txt_sce <- filterPixels(txt_sce, minevents = 40, correct_pixels = TRUE)
  txt_sce <- computeSpillmat(txt_sce)
  
  assay(cur_sce, "exprs") <- asinh(counts(cur_sce) / 5)
  bc_key_tiff <- as.numeric(unique(cur_sce$sample_mass))
  bc_key_tiff <- bc_key_tiff[order(bc_key_tiff)]
  cur_sce <- assignPrelim(cur_sce, bc_key = bc_key_tiff)
  cur_sce <- estCutoffs(cur_sce)
  cur_sce <- applyCutoffs(cur_sce)
  cur_sce <- filterPixels(cur_sce, minevents = 40, correct_pixels = TRUE)
  cur_sce <- computeSpillmat(cur_sce)
  
  expect_equal(
    metadata(cur_sce)$spillover_matrix,
    metadata(txt_sce)$spillover_matrix,
    tolerance = 1e-3
  )
  
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
