test_that("read_cpout function works.", {
    path <- system.file("extdata/mockData/cpout", package = "imcRtools")

    # SpatialExperiment
    cur_spe <- read_cpout(path)
    
    expect_s4_class(cur_spe, "SpatialExperiment")
    
    expect_equal(rownames(cur_spe), c("Ag107", "Pr141", "Sm147", 
                                      "Eu153", "Yb172"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 409))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Neighbors_NumberOfNeighbors_8", 
                                            "Metadata_acname", 
                                            "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStackFiltered", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- merge(object_file, image_file, by = "ImageNumber", order = FALSE)
    
    expect_equal(counts(cur_spe), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_spe)[1:10], c(0.03703704, 0.05555556, 0.05555556, 0.09259259, 0.62962963,
                                          0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.33333333))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_spe$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_spe$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_spe$Metadata_description, object_file$Metadata_description)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$full, cur_panel$full[order(mass, decreasing = FALSE)])
    
    expect_equal(colPairNames(cur_spe), "neighbourhood")
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighbourhood"))
    cur_test <- vroom::vroom(file.path(path, "Object relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 3,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_spe$sample_id %in% c("1", "2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "3")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_spe$sample_id %in% c("1", "2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "3")])
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighbourhood"))
    cur_test <- vroom::vroom(file.path(path, "Object relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 4,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_spe$sample_id %in% c("1", "2", "3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "4")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_spe$sample_id %in% c("1", "2", "3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "4")])
    
    # SingleCellExperiment
    cur_sce <- read_cpout(path, return_as = "sce")
    
    expect_s4_class(cur_sce, "SingleCellExperiment")
    
    expect_equal(rownames(cur_sce), c("Ag107", "Pr141", "Sm147", 
                                      "Eu153", "Yb172"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 409))
    expect_equal(names(rowData(cur_sce)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", "Pos_Y",
                                            "Neighbors_NumberOfNeighbors_8", 
                                            "Metadata_acname", 
                                            "Metadata_acid", "Metadata_description"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStackFiltered", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_sce)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- merge(object_file, image_file, by = "ImageNumber", order = FALSE)
    
    expect_equal(counts(cur_sce), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_sce)[1:10], c(0.03703704, 0.05555556, 0.05555556, 0.09259259, 0.62962963,
                                          0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.33333333))
    
    expect_equal(cur_sce$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_sce$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_sce$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_sce$Metadata_description, object_file$Metadata_description)
    expect_equal(cur_sce$Pos_X, object_file$Location_Center_X)
    expect_equal(cur_sce$Pos_Y, object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_sce)$Metal.Tag, cur_panel$`Metal Tag`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_sce)$Target, cur_panel$Target[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_sce)$full, cur_panel$full[order(mass, decreasing = FALSE)])
    
    expect_equal(colPairNames(cur_sce), "neighbourhood")
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighbourhood"))
    cur_test <- vroom::vroom(file.path(path, "Object relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 3,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_sce$sample_id %in% c("1", "2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "3")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_sce$sample_id %in% c("1", "2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "3")])
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighbourhood"))
    cur_test <- vroom::vroom(file.path(path, "Object relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 4,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_sce$sample_id %in% c("1", "2", "3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "4")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_sce$sample_id %in% c("1", "2", "3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "4")])
    
    # Test other inputs
    cur_spe <- read_cpout(path, image_file = NULL, scale_intensities = FALSE)
    
    expect_equal(rownames(cur_spe), c("Ag107", "Pr141", "Sm147", 
                                      "Eu153", "Yb172"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 409))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Neighbors_NumberOfNeighbors_8"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStackFiltered", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    expect_equal(counts(cur_spe), cur_counts)
    expect_equal(counts(cur_spe)[1:10], c(5.651490e-07, 8.477234e-07, 8.477234e-07, 1.412872e-06, 9.607532e-06,
                                          0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.086341e-06))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$full, cur_panel$full[order(mass, decreasing = FALSE)])
    
    expect_equal(colPairNames(cur_spe), "neighbourhood")
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighbourhood"))
    cur_test <- vroom::vroom(file.path(path, "Object relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 3,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_spe$sample_id %in% c("1", "2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "3")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_spe$sample_id %in% c("1", "2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "3")])
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighbourhood"))
    cur_test <- vroom::vroom(file.path(path, "Object relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 4,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_spe$sample_id %in% c("1", "2", "3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "4")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_spe$sample_id %in% c("1", "2", "3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "4")])
    
    cur_sce <- read_cpout(path, return_as = "sce", image_file = NULL, scale_intensities = FALSE)
    
    expect_equal(rownames(cur_sce), c("Ag107", "Pr141", "Sm147", 
                                      "Eu153", "Yb172"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 409))
    expect_equal(names(rowData(cur_sce)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", "Pos_Y",
                                            "Neighbors_NumberOfNeighbors_8"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStackFiltered", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_sce)

    expect_equal(counts(cur_sce), cur_counts)
    expect_equal(counts(cur_sce)[1:10], c(5.651490e-07, 8.477234e-07, 8.477234e-07, 1.412872e-06, 9.607532e-06,
                                          0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.086341e-06))
    
    expect_equal(cur_sce$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_sce$Pos_X, object_file$Location_Center_X)
    expect_equal(cur_sce$Pos_Y, object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_sce)$Metal.Tag, cur_panel$`Metal Tag`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_sce)$Target, cur_panel$Target[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_sce)$full, cur_panel$full[order(mass, decreasing = FALSE)])
    
    expect_equal(colPairNames(cur_sce), "neighbourhood")
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighbourhood"))
    cur_test <- vroom::vroom(file.path(path, "Object relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 3,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_sce$sample_id %in% c("1", "2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "3")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_sce$sample_id %in% c("1", "2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "3")])
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighbourhood"))
    cur_test <- vroom::vroom(file.path(path, "Object relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 4,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_sce$sample_id %in% c("1", "2", "3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "4")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_sce$sample_id %in% c("1", "2", "3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "4")])
    
    cur_spe <- read_cpout(path, scale_intensities = FALSE)
    
    expect_equal(rownames(cur_spe), c("Ag107", "Pr141", "Sm147", 
                                      "Eu153", "Yb172"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 409))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Neighbors_NumberOfNeighbors_8",
                                            "Metadata_acname", "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStackFiltered", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    expect_equal(counts(cur_spe), cur_counts)
    expect_equal(counts(cur_spe)[1:10], c(5.651490e-07, 8.477234e-07, 8.477234e-07, 1.412872e-06, 9.607532e-06,
                                          0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 5.086341e-06))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$full, cur_panel$full[order(mass, decreasing = FALSE)])
    
    expect_equal(colPairNames(cur_spe), "neighbourhood")
    
    cur_spe <- read_cpout(path, graph_file = NULL)
    
    expect_equal(rownames(cur_spe), c("Ag107", "Pr141", "Sm147", 
                                      "Eu153", "Yb172"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 409))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Neighbors_NumberOfNeighbors_8", 
                                            "Metadata_acname", 
                                            "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStackFiltered", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- merge(object_file, image_file, by = "ImageNumber", order = FALSE)
    
    expect_equal(counts(cur_spe), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_spe)[1:10], c(0.03703704, 0.05555556, 0.05555556, 0.09259259, 0.62962963,
                                          0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.33333333))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_spe$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_spe$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_spe$Metadata_description, object_file$Metadata_description)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$full, cur_panel$full[order(mass, decreasing = FALSE)])
    
    expect_equal(length(colPairNames(cur_spe)), 0)
    
    cur_spe <- read_cpout(path, panel_file = NULL)
    
    expect_equal(rownames(cur_spe), c("Ag107", "Pr141", "Sm147", 
                                      "Eu153", "Yb172"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 409))
    expect_equal(length(names(rowData(cur_spe))), 0) 
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Neighbors_NumberOfNeighbors_8", 
                                            "Metadata_acname", 
                                            "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStackFiltered", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- merge(object_file, image_file, by = "ImageNumber", order = FALSE)
    
    expect_equal(counts(cur_spe), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_spe)[1:10], c(0.03703704, 0.05555556, 0.05555556, 0.09259259, 0.62962963,
                                          0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.33333333))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_spe$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_spe$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_spe$Metadata_description, object_file$Metadata_description)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    expect_equal(colPairNames(cur_spe), "neighbourhood")
    
    cur_spe <- read_cpout(path, intensities = "MedianIntensity_FullStack_")

    expect_equal(rownames(cur_spe), c("Ag107", "Pr141", "Sm147", 
                                      "Eu153", "Yb172"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 409))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Neighbors_NumberOfNeighbors_8", 
                                            "Metadata_acname", 
                                            "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MedianIntensity_FullStack_", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- merge(object_file, image_file, by = "ImageNumber", order = FALSE)
    
    expect_equal(counts(cur_spe), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_spe)[1:10], c(0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0.5))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_spe$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_spe$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_spe$Metadata_description, object_file$Metadata_description)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik[order(mass, decreasing = FALSE)])
    expect_equal(rowData(cur_spe)$full, cur_panel$full[order(mass, decreasing = FALSE)])
    
    expect_equal(colPairNames(cur_spe), "neighbourhood")
    
    cur_spe <- read_cpout(path, extract_cellmetadata_from = "Location_MaxIntensity_Y_ProbSegmentation_c1")
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Location_MaxIntensity_Y_ProbSegmentation_c1", 
                                            "Metadata_acname", 
                                            "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    
    expect_equal(cur_spe$Location_MaxIntensity_Y_ProbSegmentation_c1, object_file$Location_MaxIntensity_Y_ProbSegmentation_c1)

    cur_spe <- read_cpout(path, extract_imagemetadata_from = "Metadata_end_timestamp")
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Neighbors_NumberOfNeighbors_8", 
                                            "Metadata_end_timestamp"))
    
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStackFiltered", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- merge(object_file, image_file, by = "ImageNumber", order = FALSE)
    
    expect_equal(cur_spe$Metadata_end_timestamp, object_file$Metadata_end_timestamp)
    
    cur_spe <- read_cpout(path, extract_coords_from = NULL)
    
    cur_spe <- read_cpout(path, extract_cellmetadata_from = NULL)
    
    cur_spe <- read_cpout(path, extract_imagemetadata_from = NULL)
    
    
})
