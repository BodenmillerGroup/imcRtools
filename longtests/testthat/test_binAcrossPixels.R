test_that("binAcrossPixels function works.", {
    path <- system.file("extdata/spillover", package = "imcRtools")
    
    # Read in .txt
    expect_silent(cur_sce <- readSCEfromTXT(path, verbose = FALSE))
    
    # Works
    expect_silent(out <- binAcrossPixels(cur_sce, bin_size = 10))
    expect_s4_class(out, "SingleCellExperiment")
    expect_equal(dim(out), c(nrow(cur_sce), ncol(cur_sce)/10))
    
    # colData are correct
    expect_equal(out$sample_id, rep(c("Dy161", "Dy162", "Dy163", "Dy164"), each = 10))
    expect_equal(out$spot_id, rep(c("Dy161", "Dy162", "Dy163", "Dy164"), each = 10))
    expect_equal(as.numeric(out$bin), rep(1:10, 4))
    expect_equal(as.numeric(out$ncells), rep(10, 40))
    
    # rowData are correct
    expect_equal(rowData(cur_sce), rowData(out))
    
    # Summarized counts are correct
    expect_equal(counts(out)[,1], rowSums(counts(cur_sce)[,1:10]))
    expect_equal(counts(out)[,2], rowSums(counts(cur_sce)[,11:20]))
    expect_equal(counts(out)[,out$sample_id == "Dy164"][,1], rowSums(counts(cur_sce)[,cur_sce$sample_id == "Dy164"][,1:10]))
    expect_equal(counts(out)[,out$sample_id == "Dy164"][,2], rowSums(counts(cur_sce)[,cur_sce$sample_id == "Dy164"][,11:20]))
    
    # Works
    expect_silent(out <- binAcrossPixels(cur_sce, bin_size = 2))
    expect_s4_class(out, "SingleCellExperiment")
    expect_equal(dim(out), c(nrow(cur_sce), ncol(cur_sce)/2))
    
    # colData are correct
    expect_equal(out$sample_id, rep(c("Dy161", "Dy162", "Dy163", "Dy164"), each = 50))
    expect_equal(out$spot_id, rep(c("Dy161", "Dy162", "Dy163", "Dy164"), each = 50))
    expect_equal(as.numeric(out$bin), rep(1:50, 4))
    expect_equal(as.numeric(out$ncells), rep(2, 200))
    
    # rowData are correct
    expect_equal(rowData(cur_sce), rowData(out))
    
    # Summarized counts are correct
    expect_equal(counts(out)[,1], rowSums(counts(cur_sce)[,1:2]))
    expect_equal(counts(out)[,2], rowSums(counts(cur_sce)[,3:4]))
    expect_equal(counts(out)[,3], rowSums(counts(cur_sce)[,5:6]))
    expect_equal(counts(out)[,4], rowSums(counts(cur_sce)[,7:8]))
    expect_equal(counts(out)[,out$sample_id == "Dy164"][,1], rowSums(counts(cur_sce)[,cur_sce$sample_id == "Dy164"][,1:2]))
    expect_equal(counts(out)[,out$sample_id == "Dy164"][,2], rowSums(counts(cur_sce)[,cur_sce$sample_id == "Dy164"][,3:4]))
    expect_equal(counts(out)[,out$sample_id == "Dy164"][,3], rowSums(counts(cur_sce)[,cur_sce$sample_id == "Dy164"][,5:6]))
    expect_equal(counts(out)[,out$sample_id == "Dy164"][,4], rowSums(counts(cur_sce)[,cur_sce$sample_id == "Dy164"][,7:8]))
    
    cur_sce$test_id <- cur_sce$sample_id
    expect_silent(out <- binAcrossPixels(cur_sce, bin_size = 5, spot_id = "test_id"))
    
    assay(cur_sce, "test") <- counts(cur_sce)
    expect_silent(out <- binAcrossPixels(cur_sce, bin_size = 5, assay_type = "test"))
    
    expect_silent(out <- binAcrossPixels(cur_sce, bin_size = 5, statistic = "mean"))
    
    expect_equal(counts(out)[,1], rowMeans(counts(cur_sce)[,1:5]))
    expect_equal(counts(out)[,2], rowMeans(counts(cur_sce)[,6:10]))
    expect_equal(counts(out)[,out$sample_id == "Dy164"][,1], rowMeans(counts(cur_sce)[,cur_sce$sample_id == "Dy164"][,1:5]))
    expect_equal(counts(out)[,out$sample_id == "Dy164"][,2], rowMeans(counts(cur_sce)[,cur_sce$sample_id == "Dy164"][,6:10]))
    
    # Error
    expect_error(binAcrossPixels("test"),
                 regexp = "'object' needs to be a SingleCellExperiment object.",
                 fixed = TRUE)
    expect_error(binAcrossPixels(cur_sce, bin_size = c(1,2)),
                 regexp = "'bin_size' needs to be single numeric.",
                 fixed = TRUE)
    expect_error(binAcrossPixels(cur_sce, bin_size = "test"),
                 regexp = "'bin_size' needs to be single numeric.",
                 fixed = TRUE)
    expect_error(binAcrossPixels(cur_sce, spot_id = "test"),
                 regexp = "'spot_id' not in 'colData(object)'.",
                 fixed = TRUE)
    expect_error(binAcrossPixels(cur_sce, bin_size = 10, assay_type = "test_1"),
                 regexp = "'assay_type' not in 'assayNames(object)'.",
                 fixed = TRUE)
    expect_error(binAcrossPixels(cur_sce, bin_size = 10, statistic = "test"),
                 regexp = "'statistic' must be 'sum', 'mean' or 'median'.",
                 fixed = TRUE)
})
