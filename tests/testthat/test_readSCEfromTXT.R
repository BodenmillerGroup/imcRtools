test_that("readSCEfromTXT function reads in correct objects.", {
    path <- system.file("extdata/spillover", package = "imcRtools")
    
    cur_files <- list.files(path, pattern = ".txt", full.names = TRUE)
    
    # Read in with read.csv
    cur_files_txt <- lapply(cur_files, read.delim)
    names(cur_files_txt)
})
