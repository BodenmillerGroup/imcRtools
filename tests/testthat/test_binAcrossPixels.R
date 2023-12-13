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
    
    test <- aggregate(t(counts(cur_sce)), by = list(cur_sce$sample_id, rep(rep(1:10, each = 10), 4)), FUN = "sum")
    test <- test[order(test$Group.1, test$Group.2),]
    
    test2 <- counts(out)
    test <- t(as.matrix(test[,-c(1,2)]))
    
    dimnames(test2) <- NULL
    dimnames(test)  <- NULL
    
    expect_equal(test2, test)
    
    expect_silent(out <- binAcrossPixels(cur_sce, bin_size = 10, statistic = "mean"))
    
    test <- aggregate(t(counts(cur_sce)), by = list(cur_sce$sample_id, rep(rep(1:10, each = 10), 4)), FUN = "mean")
    test <- test[order(test$Group.1, test$Group.2),]
    
    test2 <- counts(out)
    test <- t(as.matrix(test[,-c(1,2)]))
    
    dimnames(test2) <- NULL
    dimnames(test)  <- NULL
    
    expect_equal(test2, test)
    
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
    
    # Works
    expect_silent(out <- binAcrossPixels(cur_sce, bin_size = 3))
    expect_s4_class(out, "SingleCellExperiment")
    
    # colData are correct
    expect_equal(out$sample_id, rep(c("Dy161", "Dy162", "Dy163", "Dy164"), each = 34))
    expect_equal(out$spot_id, rep(c("Dy161", "Dy162", "Dy163", "Dy164"), each = 34))
    expect_equal(as.numeric(out$bin), rep(1:34, 4))
    expect_equal(as.numeric(out$ncells), rep(c(rep(3, 33), 1), 4))
    
    # rowData are correct
    expect_equal(rowData(cur_sce), rowData(out))
    
    # Summarized counts are correct
    expect_equal(counts(out)[,1], rowSums(counts(cur_sce)[,1:3]))
    expect_equal(counts(out)[,2], rowSums(counts(cur_sce)[,4:6]))
    expect_equal(counts(out)[,3], rowSums(counts(cur_sce)[,7:9]))
    expect_equal(counts(out)[,4], rowSums(counts(cur_sce)[,10:12]))
    expect_equal(counts(out)[,34], counts(cur_sce)[,100])    
    
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
    
    cur_sce <- cur_sce[,sample(ncol(cur_sce))]
    
    expect_error(binAcrossPixels(cur_sce, bin_size = 10),
                 regexp = "Spot IDs of pixels within 'object' are not ordered alphabetically.",
                 fixed = TRUE)
})
