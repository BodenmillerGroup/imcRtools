test_that("readImagefromTXT function works.", {
    path <- system.file("extdata/mockData/raw", package = "imcRtools")
    
    # Works 
    expect_silent(cur_cil <- readImagefromTXT(path))
    
    expect_s4_class(cur_cil, "CytoImageList")
    expect_equal(length(cur_cil), 3)
    expect_equal(channelNames(cur_cil), c("Ag107Di", "Pr141Di", "Sm147Di", "Eu153Di", "Yb172Di"))
    expect_equal(names(cur_cil), c("20210305_NE_mockData2_ROI_001_1",  
                                   "20210305_NE_mockData2_ROI_002_2", 
                                   "20210305_NE_mockData2_ROI_003_3"))
    
    # Test for numerical correctness
    cur_tiff <- system.file("extdata/mockData/tiffs/20210305_NE_mockData2_s0_a1_ac_fullFiltered.tiff", package = "imcRtools")
    test_cil <- cytomapper::loadImages(cur_tiff, as.is = TRUE)
    
    expect_equal(dim(cur_cil[[1]]), dim(test_cil[[1]]))
    test_1 <- floor(as.array(cur_cil[[1]][,,1]))
    dimnames(test_1) <- NULL
    test_2 <- as.array(test_cil[[1]][,,1])
    expect_equal(test_1,
                 test_2)
    
    test_1 <- floor(as.array(cur_cil[[1]][,,2]))
    dimnames(test_1) <- NULL
    test_2 <- as.array(test_cil[[1]][,,2])
    expect_equal(test_1,
                 test_2)
    
    test_1 <- floor(as.array(cur_cil[[1]][,,2]))
    dimnames(test_1) <- NULL
    test_2 <- as.array(test_cil[[1]][,,2])
    expect_equal(test_1,
                 test_2)
    
    test_1 <- floor(as.array(cur_cil[[1]][,,2]))
    dimnames(test_1) <- NULL
    test_2 <- as.array(test_cil[[1]][,,2])
    expect_equal(test_1,
                 test_2)
    
    test_1 <- floor(as.array(cur_cil[[1]][,,2]))
    dimnames(test_1) <- NULL
    test_2 <- as.array(test_cil[[1]][,,2])
    expect_equal(test_1,
                 test_2)
    
    # Read in individual files
    expect_silent(cur_cil <- readImagefromTXT(path, pattern = "ROI_002"))
    expect_equal(length(cur_cil), 1)
    expect_equal(channelNames(cur_cil), c("Ag107Di", "Pr141Di", "Sm147Di", "Eu153Di", "Yb172Di"))
    expect_equal(names(cur_cil), c("20210305_NE_mockData2_ROI_002_2"))
    
    # Read in different channelNames
    expect_silent(cur_cil <- readImagefromTXT(path, channel_pattern = "[A-Za-z]{2}[0-9]{3}"))
    expect_equal(length(cur_cil), 3)
    expect_equal(channelNames(cur_cil), c("Ag107", "Pr141", "Sm147", "Eu153", "Yb172"))
    expect_equal(names(cur_cil), c("20210305_NE_mockData2_ROI_001_1",  
                                   "20210305_NE_mockData2_ROI_002_2", 
                                   "20210305_NE_mockData2_ROI_003_3"))
    
    # Read in single channel
    expect_silent(cur_cil <- readImagefromTXT(path, channel_pattern = "Ag107"))
    expect_equal(length(cur_cil), 3)
    expect_equal(channelNames(cur_cil), c("Ag107"))
    expect_equal(names(cur_cil), c("20210305_NE_mockData2_ROI_001_1",  
                                   "20210305_NE_mockData2_ROI_002_2", 
                                   "20210305_NE_mockData2_ROI_003_3"))
    
    test_1 <- floor(as.array(cur_cil[[1]][,,1]))
    dimnames(test_1) <- NULL
    test_2 <- as.array(test_cil[[1]][,,1])
    expect_equal(test_1,
                 test_2)
    
    # parallelisation
    expect_silent(cur_cil <- readImagefromTXT(path, 
                                              BPPARAM = BiocParallel::bpparam()))
    expect_equal(length(cur_cil), 3)
    expect_equal(channelNames(cur_cil), c("Ag107Di", "Pr141Di", "Sm147Di", "Eu153Di", "Yb172Di"))
    expect_equal(names(cur_cil), c("20210305_NE_mockData2_ROI_001_1",  
                                   "20210305_NE_mockData2_ROI_002_2", 
                                   "20210305_NE_mockData2_ROI_003_3"))
    
    # Error
    expect_error(cur_cil <- readImagefromTXT("test"),
                 regexp = "The indicated path does not exist",
                 fixed = TRUE)
    
    expect_error(cur_cil <- readImagefromTXT(1),
                 regexp = "Please provide a string input indicating a path.",
                 fixed = TRUE)
    
    expect_error(cur_cil <- readImagefromTXT(path, pattern = 1),
                 regexp = "'pattern' must belibrar indicated as single character",
                 fixed = TRUE)
    
    expect_error(cur_cil <- readImagefromTXT(path, pattern = "test"),
                 regexp = "The pattern does not match any\nof the files in the provided directory.",
                 fixed = TRUE)
    
    expect_error(cur_cil <- readImagefromTXT(path, index_names = c("test", "test2")),
                 regexp = "'index_names' not in the names of the .txt files.",
                 fixed = TRUE)
    
    expect_error(cur_cil <- readImagefromTXT(path, channel_pattern = "test"),
                 regexp = "'channel_pattern' does not match any entries in the .txt files.",
                 fixed = TRUE)
    
})
