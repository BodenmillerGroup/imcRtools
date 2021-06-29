test_that("read_steinbock and read_cpout generate the same objects", {
    steinbock_path <- "Downloads/with_mean_factor100/"
    cpout_path <- system.file("extdata/mockData/cpout", package = "imcRtools")
  
    # SpatialExperiment
    steinbock_spe <- read_steinbock(steinbock_path)
    cpout_spe <- read_cpout(cpout_path)
    
    
})
