test_that("findBorderCells function works", {
    path <- system.file("extdata/mockData/steinbock", package = "imcRtools")
    spe <- read_steinbock(path)

    spe  <- findBorderCells(spe, img_id = "sample_id", border_dist = 10)
    
    plotSpatial(spe, img_id = "sample_id", node_color_by = "border_cells")
    
    cur_ind <- unique(spe$sample_id)
    
    out <- lapply(cur_ind, function(x){
        cur_obj <- spe[,spe$sample_id == x]
        min_x <- min(spatialCoords(cur_obj)[,"Pos_X"])
        max_x <- max(spatialCoords(cur_obj)[,"Pos_X"])
        min_y <- min(spatialCoords(cur_obj)[,"Pos_Y"])
        max_y <- max(spatialCoords(cur_obj)[,"Pos_Y"])
        cur_x <- spatialCoords(cur_obj)[,"Pos_X"]
        cur_y <- spatialCoords(cur_obj)[,"Pos_Y"]
        cur_cells <- cur_x <= min_x + 10 | cur_x >= max_x - 10 |
            cur_y <= min_y + 10 | cur_y >= max_y - 10
        cur_obj$border_cells <- cur_cells
        return(cur_obj)
    })
    
    out <- do.call("cbind", out)
    
    expect_equal(spe$border_cells, out$border_cells, check.attributes = FALSE)
    
    spe  <- findBorderCells(spe, img_id = "sample_id", border_dist = 5)
    
    plotSpatial(spe, img_id = "sample_id", node_color_by = "border_cells")
    
    cur_ind <- unique(spe$sample_id)
    
    out <- lapply(cur_ind, function(x){
        cur_obj <- spe[,spe$sample_id == x]
        min_x <- min(spatialCoords(cur_obj)[,"Pos_X"])
        max_x <- max(spatialCoords(cur_obj)[,"Pos_X"])
        min_y <- min(spatialCoords(cur_obj)[,"Pos_Y"])
        max_y <- max(spatialCoords(cur_obj)[,"Pos_Y"])
        cur_x <- spatialCoords(cur_obj)[,"Pos_X"]
        cur_y <- spatialCoords(cur_obj)[,"Pos_Y"]
        cur_cells <- cur_x <= min_x + 5 | cur_x >= max_x - 5 |
            cur_y <= min_y + 5 | cur_y >= max_y - 5
        cur_obj$border_cells <- cur_cells
        return(cur_obj)
    })
    
    out <- do.call("cbind", out)
    
    expect_equal(spe$border_cells, out$border_cells, check.attributes = FALSE)
    
    path <- system.file("extdata/mockData/steinbock", package = "imcRtools")
    sce <- read_steinbock(path, return_as = "sce")
    
    sce  <- findBorderCells(sce, img_id = "sample_id", border_dist = 10)
    
    plotSpatial(sce, img_id = "sample_id", node_color_by = "border_cells")
    
    cur_ind <- unique(sce$sample_id)
    
    out <- lapply(cur_ind, function(x){
        cur_obj <- sce[,sce$sample_id == x]
        min_x <- min(colData(cur_obj)[,"Pos_X"])
        max_x <- max(colData(cur_obj)[,"Pos_X"])
        min_y <- min(colData(cur_obj)[,"Pos_Y"])
        max_y <- max(colData(cur_obj)[,"Pos_Y"])
        cur_x <- colData(cur_obj)[,"Pos_X"]
        cur_y <- colData(cur_obj)[,"Pos_Y"]
        cur_cells <- cur_x <= min_x + 10 | cur_x >= max_x - 10 |
            cur_y <= min_y + 10 | cur_y >= max_y - 10
        cur_obj$border_cells <- cur_cells
        return(cur_obj)
    })
    
    out <- do.call("cbind", out)
    
    expect_equal(sce$border_cells, out$border_cells, check.attributes = FALSE)
    
    sce  <- findBorderCells(sce, img_id = "sample_id", border_dist = 5)
    
    plotSpatial(sce, img_id = "sample_id", node_color_by = "border_cells")
    
    cur_ind <- unique(sce$sample_id)
    
    out <- lapply(cur_ind, function(x){
        cur_obj <- sce[,sce$sample_id == x]
        min_x <- min(colData(cur_obj)[,"Pos_X"])
        max_x <- max(colData(cur_obj)[,"Pos_X"])
        min_y <- min(colData(cur_obj)[,"Pos_Y"])
        max_y <- max(colData(cur_obj)[,"Pos_Y"])
        cur_x <- colData(cur_obj)[,"Pos_X"]
        cur_y <- colData(cur_obj)[,"Pos_Y"]
        cur_cells <- cur_x <= min_x + 5 | cur_x >= max_x - 5 |
            cur_y <= min_y + 5 | cur_y >= max_y - 5
        cur_obj$border_cells <- cur_cells
        return(cur_obj)
    })
    
    out <- do.call("cbind", out)
    
    expect_equal(sce$border_cells, out$border_cells, check.attributes = FALSE)
    
    # SingleCellExperiment
    library(cytomapper)
    data(pancreasSCE)
    
    cur_sce1 <- pancreasSCE[,pancreasSCE$ImageNb == 1]
    cur_sce2 <- pancreasSCE[,pancreasSCE$ImageNb == 2]
    cur_sce3 <- pancreasSCE[,pancreasSCE$ImageNb == 3]
    
    cur_sce1$Pos_X <- cur_sce1$Pos_X - min(cur_sce1$Pos_X)
    cur_sce1$Pos_Y <- cur_sce1$Pos_Y - min(cur_sce1$Pos_Y)
    cur_sce2$Pos_X <- cur_sce2$Pos_X - min(cur_sce2$Pos_X)
    cur_sce2$Pos_Y <- cur_sce2$Pos_Y - min(cur_sce2$Pos_Y)
    cur_sce3$Pos_X <- cur_sce3$Pos_X - min(cur_sce3$Pos_X)
    cur_sce3$Pos_Y <- cur_sce3$Pos_Y - min(cur_sce3$Pos_Y)
    
    pancreasSCE <- cbind(cur_sce1, cur_sce2, cur_sce3)

    expect_silent(pancreasSCE  <- findBorderCells(pancreasSCE, img_id = "ImageNb", border_dist = 10))
    plotSpatial(pancreasSCE, img_id = "ImageNb", node_color_by = "border_cells")
    
    # Error
    expect_error(findBorderCells("test"), 
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(findBorderCells(sce, img_id = c(1,2)), 
                 regexp = "'img_id' must be a single string.",
                 fixed = TRUE)
    expect_error(findBorderCells(sce, img_id = "test"), 
                 regexp = "'img_id' not in colData(object).",
                 fixed = TRUE)
    expect_error(findBorderCells(sce, img_id = "sample_id", coords = "test"), 
                 regexp = "'coords' must be a character vector of length 2.",
                 fixed = TRUE)
    expect_error(findBorderCells(sce, img_id = "sample_id", coords = c("Pos_x", "Pos_Y")), 
                 regexp = "'coords' not in colData(object).",
                 fixed = TRUE)
    expect_error(findBorderCells(spe, img_id = "sample_id", coords = c("Pos_x", "Pos_Y")), 
                 regexp = "'coords' not in spatialCoords(object).",
                 fixed = TRUE)
    expect_error(findBorderCells(spe, img_id = "sample_id", border_dist = "test"), 
                 regexp = "'border_dist' must be a single numeric.",
                 fixed = TRUE)
    expect_error(findBorderCells(spe, img_id = "sample_id", border_dist = c(1,2)), 
                 regexp = "'border_dist' must be a single numeric.",
                 fixed = TRUE)
})

test_that("findBorderCells function works if cells are not ordered by image", {
    # SingleCellExperiment
    library(cytomapper)
    data(pancreasSCE)
    
    cur_sce1 <- pancreasSCE[,pancreasSCE$ImageNb == 1]
    cur_sce2 <- pancreasSCE[,pancreasSCE$ImageNb == 2]
    cur_sce3 <- pancreasSCE[,pancreasSCE$ImageNb == 3]
    
    cur_sce1$Pos_X <- cur_sce1$Pos_X - min(cur_sce1$Pos_X)
    cur_sce1$Pos_Y <- cur_sce1$Pos_Y - min(cur_sce1$Pos_Y)
    cur_sce2$Pos_X <- cur_sce2$Pos_X - min(cur_sce2$Pos_X)
    cur_sce2$Pos_Y <- cur_sce2$Pos_Y - min(cur_sce2$Pos_Y)
    cur_sce3$Pos_X <- cur_sce3$Pos_X - min(cur_sce3$Pos_X)
    cur_sce3$Pos_Y <- cur_sce3$Pos_Y - min(cur_sce3$Pos_Y)
    
    pancreasSCE <- cbind(cur_sce1, cur_sce2, cur_sce3)
    
    expect_silent(pancreasSCE  <- findBorderCells(pancreasSCE, img_id = "ImageNb", border_dist = 10))
    plotSpatial(pancreasSCE, img_id = "ImageNb", node_color_by = "border_cells")
    
    set.seed(123)
    sce2 <- pancreasSCE[,sample(ncol(pancreasSCE))]
    
    expect_silent(sce2  <- findBorderCells(sce2, img_id = "ImageNb", border_dist = 10))
    plotSpatial(sce2, img_id = "ImageNb", node_color_by = "border_cells")
    
    expect_equal(pancreasSCE$border_cells, sce2[,colnames(pancreasSCE)]$border_cells)

    })
