test_that("readSCEfromTXT function reads in correct objects.", {
    path <- system.file("extdata/spillover", package = "imcRtools")
    
    # Read in .txt
    expect_silent(cur_sce <- readSCEfromTXT(path, verbose = FALSE))
    expect_equal(rowData(cur_sce)$channel_name, c("Dy161Di", "Dy162Di", 
                                                  "Dy163Di","Dy164Di"))
    expect_equal(rowData(cur_sce)$marker_name, c("Dy161", "Dy162", 
                                                  "Dy163","Dy164"))
    expect_equal(cur_sce$Start_push[1:10], c(1570, 1955, 2341, 2726, 3111,
                                             3497, 3882, 4267, 4653, 5038))
    expect_equal(cur_sce$End_push[1:10], c(1953, 2338, 2724, 3109, 3494, 
                                           3880, 4265, 4650, 5036, 5421))
    expect_equal(cur_sce$Pushes_duration[1:10], cur_sce$End_push[1:10] - cur_sce$Start_push[1:10] + 1)
    expect_equal(cur_sce$X, rep(0:99, 4))
    expect_equal(cur_sce$Y, rep(0, 400))
    expect_equal(cur_sce$Z, cur_sce$X)
    expect_equal(cur_sce$sample_id, rep(c("Dy161", "Dy162", 
                                          "Dy163","Dy164"), each = 100))
    expect_equal(cur_sce$sample_metal, rep(c("Dy", "Dy", 
                                          "Dy","Dy"), each = 100))
    expect_equal(cur_sce$sample_mass, rep(c("161", "162", 
                                            "163","164"), each = 100))
    
    expect_equal(assayNames(cur_sce), "counts")
    
    cur_txt <- lapply(list.files(path, full.names = TRUE), read_delim, delim = "\t")
    cur_txt <- lapply(cur_txt, as.data.frame)
    cur_txt <- do.call("rbind", cur_txt)
    
    expect_equal(as.numeric(counts(cur_sce)[1,]), cur_txt$`161Dy(Dy161Di)`)
    expect_equal(as.numeric(counts(cur_sce)[2,]), cur_txt$`162Dy(Dy162Di)`)
    expect_equal(as.numeric(counts(cur_sce)[3,]), cur_txt$`163Dy(Dy163Di)`)
    expect_equal(as.numeric(counts(cur_sce)[4,]), cur_txt$`164Dy(Dy164Di)`)
    expect_equal(dim(cur_sce), c(4, 400))
    
    # Verbose output
    cur_out <- capture_output(cur_sce <- readSCEfromTXT(path))
    expect_equal(cur_out, "Spotted channels:  Dy161, Dy162, Dy163, Dy164\nAcquired channels:  Dy161, Dy162, Dy163, Dy164\nChannels spotted but not acquired:  \nChannels acquired but not spotted:  ")
    
    # Other parameters
    expect_silent(cur_sce_2 <- readSCEfromTXT(path, pattern = "Dy162", verbose = FALSE))
    expect_equal(dim(cur_sce_2), c(4, 100))
    expect_equal(rowData(cur_sce)$channel_name, c("Dy161Di", "Dy162Di", 
                                                  "Dy163Di","Dy164Di"))
    
    expect_silent(cur_sce_2 <- readSCEfromTXT(path, metadata_cols = "X", verbose = FALSE))
    expect_equal(counts(cur_sce), counts(cur_sce_2))
    expect_equal(names(colData(cur_sce_2)), c("X", "sample_id", "sample_metal", 
                                              "sample_mass" ))
    
    # Read in list
    cur_files <- list.files(path, pattern = ".txt", full.names = TRUE)
    cur_files_names <- list.files(path, pattern = ".txt")
    cur_files <- lapply(cur_files, read_delim, delim = "\t")
    names(cur_files) <- str_extract(cur_files_names, "[A-Za-z]{1,2}[0-9]{2,3}")
    
    expect_silent(cur_sce_3 <- readSCEfromTXT(cur_files, verbose = FALSE))
    expect_equal(cur_sce, cur_sce_3)
    expect_silent(cur_sce_4 <- readSCEfromTXT(cur_files, metadata_cols = "X", verbose = FALSE))
    expect_equal(cur_sce_2, cur_sce_4)
    
    # Don't extract metal names
    cur_sce_5 <- readSCEfromTXT(path, verbose = FALSE, read_metal_from_filename = FALSE)
    
    expect_equal(counts(cur_sce), counts(cur_sce_5))
    expect_equal(names(colData(cur_sce_5)), c("Start_push", "End_push",
                                              "Pushes_duration", "X", "Y",
                                              "Z", "sample_id"))

    # Error
    names(cur_files)[1] <- "test" 
    expect_error(readSCEfromTXT(cur_files, verbose = FALSE), 
                 regexp = "Not all names match the pattern (mt)(mass).",
                 fixed = TRUE)
    
    names(cur_files)[1] <- "XYZ1" 
    expect_error(readSCEfromTXT(cur_files, verbose = FALSE), 
                 regexp = "Not all names match the pattern (mt)(mass).",
                 fixed = TRUE)
    
    names(cur_files) <- NULL
    expect_error(readSCEfromTXT(cur_files, verbose = FALSE), 
                 regexp = "If 'x' is a list, it needs to be named.",
                 fixed = TRUE)
    
    expect_error(readSCEfromTXT("test", verbose = FALSE), 
                 regexp = "Path does not exist.",
                 fixed = TRUE)
    
    expect_error(readSCEfromTXT(c("test", "test2"), verbose = FALSE), 
                 regexp = "Input 'x' is not of the correct format.",
                 fixed = TRUE)
    
    expect_error(readSCEfromTXT(path, metadata_cols = "test", verbose = FALSE), 
                 regexp = "Not all 'metadata_cols' are present in entries to 'txt_list'",
                 fixed = TRUE)
    
})
