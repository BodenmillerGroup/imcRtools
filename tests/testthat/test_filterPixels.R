test_that("filterPixels function works.", {
    path <- system.file("extdata/spillover", package = "imcRtools")
    
    # Read in .txt
    expect_silent(cur_sce <- readSCEfromTXT(path, verbose = FALSE))
    assay(cur_sce, "exprs") <- asinh(counts(cur_sce)/5)

    bc_key <- as.numeric(unique(cur_sce$sample_mass))
    
    expect_silent(cur_sce <- assignPrelim(cur_sce, bc_key = bc_key, 
                                          verbose = FALSE))
    
    expect_true(all(is.character(cur_sce$bc_id)))
    
    expect_equal(cur_sce$sample_mass, cur_sce$bc_id)
    expect_equal(cur_sce$delta[1:5], 
                 c(0.5496888, 0.5281406, 0.4944301, 0.4765845, 0.4958716),
                 tolerance = 0.000001)
    
    cur_sce <- estCutoffs(cur_sce)
    
    expect_equal(metadata(cur_sce)$sep_cutoffs, 
                 c("161" = 0.3900550, "162" = 0.2900170, "163" = 0.3900138, "164" = 0.3500177),
                 tolerance = 0.000001)
    
    expect_silent(cur_sce <- applyCutoffs(cur_sce))
    expect_equal(sum(cur_sce$bc_id == "0"), 2)
    
    expect_true(all(is.character(cur_sce$bc_id)))
    
    # Filter function shouldn't do anything
    expect_silent(cur_sce_2 <- filterPixels(cur_sce))
    expect_identical(cur_sce_2, cur_sce)
    
    sm_1 <- CATALYST::computeSpillmat(cur_sce)
    sm_1 <- metadata(sm_1)$spillover_matrix
    
    cur_sm <- matrix(c(1.000000000, 0.031443129, 0.009734712, 0.006518048,
                       0.015715159, 1.000000000, 0.048116187, 0.008250039,
                       0.003809504, 0.012159704, 1.000000000, 0.020214651,
                       0.005058069, 0.008457546, 0.028912343, 1.000000000),
                     ncol = 4, nrow = 4, byrow = TRUE, 
                     dimnames = list(c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di"),
                                     c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")))
    
    expect_equal(sm_1, cur_sm)
    
    # Change a few pixels
    cur_sce_test <- cur_sce
    
    set.seed(12345)
    cur_sam <- sample(1:ncol(cur_sce), 20)
    
    cur_sce_test$bc_id[cur_sam] <- "test"
    expect_silent(cur_sce_2 <- filterPixels(cur_sce_test))
    
    expect_true(all(cur_sce_2$bc_id[cur_sam] == "0"))
    expect_equal(sum(cur_sce_2$bc_id == "0"), 22)
    
    cur_sce_test <- cur_sce
    cur_sce_test$sample_mass[cur_sam] <- "test"
    expect_silent(cur_sce_2 <- filterPixels(cur_sce_test))
    
    expect_true(all(cur_sce_2$bc_id[cur_sam] == "0"))
    expect_equal(sum(cur_sce_2$bc_id == "0"), 22)
    
    sm_2 <- CATALYST::computeSpillmat(cur_sce_2)
    sm_2 <- metadata(sm_2)$spillover_matrix
    
    cur_sm <- matrix(c(1.000000000, 0.031379344, 0.009734712, 0.006434588,
                       0.015716916, 1.000000000, 0.048159650, 0.008242979,
                       0.003791071, 0.012178834, 1.000000000, 0.020185447,
                       0.005052712, 0.008413039, 0.028952340, 1.000000000),
                     ncol = 4, nrow = 4, byrow = TRUE, 
                     dimnames = list(c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di"),
                                     c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")))
    
    expect_equal(sm_2, cur_sm)
    
    # Don't correct pixels
    cur_sce_test <- cur_sce
    
    set.seed(12345)
    cur_sam <- sample(1:ncol(cur_sce), 20)
    
    cur_sce_test$bc_id[cur_sam] <- "test"
    expect_silent(cur_sce_2 <- filterPixels(cur_sce_test, 
                                            correct_pixels = FALSE,
                                            minevents = 18))
    
    expect_true(all(cur_sce_2$bc_id[cur_sam] == "test"))
    
    # Filter based on size
    expect_silent(cur_sce_2 <- filterPixels(cur_sce_test, 
                                            correct_pixels = FALSE,
                                            minevents = 22))
    
    expect_true(all(cur_sce_2$bc_id[cur_sam] == "0"))
    
    cur_sce_test <- cur_sce
    expect_silent(cur_sce_2 <- filterPixels(cur_sce_test, minevents = 120))
    
    expect_true(all(cur_sce_2$bc_id == "0"))
    
    sm_3 <- CATALYST::computeSpillmat(cur_sce_2)
    sm_3 <- metadata(sm_3)$spillover_matrix
    
    cur_sm <- matrix(c(1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1),
                     ncol = 4, nrow = 4, byrow = TRUE, 
                     dimnames = list(c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di"),
                                     c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")))
    
    expect_equal(sm_3, cur_sm)
    
    cur_sce_test <- cur_sce
    cur_sce_test$bc_id[cur_sce_test$bc_id == "161"] <- c(rep("161", 20), rep("0", 80))
    expect_silent(cur_sce_2 <- filterPixels(cur_sce_test, minevents = 21))
    
    expect_true(all(cur_sce_2$bc_id[cur_sce$bc_id == "161"] == "0"))
    
    sm_4 <- CATALYST::computeSpillmat(cur_sce_2)
    sm_4 <- metadata(sm_4)$spillover_matrix
    
    cur_sm <- matrix(c(1, 0, 0, 0,
                       0.015715159, 1, 0.04811619, 0.008250039,
                       0.003809504, 0.012159704, 1, 0.020214651,
                       0.005058069, 0.008457546, 0.02891234, 1),
                     ncol = 4, nrow = 4, byrow = TRUE, 
                     dimnames = list(c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di"),
                                     c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")))
    
    expect_equal(sm_4, cur_sm)
    
    # Error
    expect_error(filterPixels("test"),
                 regexp = "'object' needs to be a SingleCellExperiment object.",
                 fixed = TRUE)
    
    expect_error(filterPixels(object = cur_sce, bc_id = "test"),
                 regexp = "'bc_id' not in 'colData(object)'",
                 fixed = TRUE)

    expect_error(filterPixels(object = cur_sce, spot_mass = "test"),
                 regexp = "'spot_mass' not in 'colData(object)'",
                 fixed = TRUE)
    
    expect_error(filterPixels(object = cur_sce, minevents = "test"),
                 regexp = "'minevents' needs to be single numeric.",
                 fixed = TRUE)
    
    expect_error(filterPixels(object = cur_sce, correct_pixels = "test"),
                 regexp = "'correct_pixels' needs to be single logical.",
                 fixed = TRUE)
})
