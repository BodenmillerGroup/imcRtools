test_that("read_cpout function works.", {

    path <- system.file("extdata/mockData/cpout", package = "imcRtools")

    # SpatialExperiment
    cur_spe <- read_cpout(path, graph_file = "Object_relationships.csv")
    
    expect_s4_class(cur_spe, "SpatialExperiment")
    
    expect_equal(rownames(cur_spe), c("Sm147", "Yb172", "Pr141", "Eu153", "Ag107"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 209))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius", "ImageNumber", "Metadata_acname", "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStack", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- dplyr::left_join(object_file, image_file, by = "ImageNumber")
    
    expect_equal(counts(cur_spe), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_spe)[1:10], c(0, 0.930232557922948, 0.0232558139480737, 0.0465116278961474, 
                                          0.116279069740369, 0.0392156862653792, 0.813725490006618, 0.107843137229793, 
                                          0.0588235293980688, 0.117647058796138))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_spe$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_spe$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_spe$Metadata_description, object_file$Metadata_description)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`)
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$full, cur_panel$full)
    
    expect_equal(colPairNames(cur_spe), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- vroom::vroom(file.path(path, "Object_relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 3,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_spe$sample_id %in% c("1", "2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "3")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_spe$sample_id %in% c("1", "2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "3")])
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- vroom::vroom(file.path(path, "Object_relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 4,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_spe$sample_id %in% c("1", "2", "3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "4")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_spe$sample_id %in% c("1", "2", "3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "4")])
    
    # SingleCellExperiment
    cur_sce <- read_cpout(path, return_as = "sce", graph_file = "Object_relationships.csv")
    
    expect_s4_class(cur_sce, "SingleCellExperiment")
    
    expect_equal(rownames(cur_sce), c("Sm147", "Yb172", "Pr141", "Eu153", "Ag107"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 209))
    expect_equal(names(rowData(cur_sce)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", "Pos_Y", "AreaShape_Area", 
                                            "Neighbors_NumberOfNeighbors_8", "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", 
                                            "AreaShape_MinorAxisLength", "AreaShape_MeanRadius", "ImageNumber", "Metadata_acname", 
                                            "Metadata_acid", "Metadata_description"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStack", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_sce)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- dplyr::left_join(object_file, image_file, by = "ImageNumber")
    
    expect_equal(counts(cur_sce), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_sce)[1:10], c(0, 0.930232557922948, 0.0232558139480737, 0.0465116278961474, 
                                          0.116279069740369, 0.0392156862653792, 0.813725490006618, 0.107843137229793, 
                                          0.0588235293980688, 0.117647058796138))
    
    expect_equal(cur_sce$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_sce$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_sce$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_sce$Metadata_description, object_file$Metadata_description)
    expect_equal(cur_sce$Pos_X, object_file$Location_Center_X)
    expect_equal(cur_sce$Pos_Y, object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`)
    expect_equal(rowData(cur_sce)$Metal.Tag, cur_panel$`Metal Tag`)
    expect_equal(rowData(cur_sce)$Target, cur_panel$Target)
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_sce)$full, cur_panel$full)
    
    expect_equal(colPairNames(cur_sce), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- vroom::vroom(file.path(path, "Object_relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 3,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_sce$sample_id %in% c("1", "2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "3")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_sce$sample_id %in% c("1", "2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "3")])
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- vroom::vroom(file.path(path, "Object_relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 4,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_sce$sample_id %in% c("1", "2", "3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "4")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_sce$sample_id %in% c("1", "2", "3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "4")])
    
    # Test other inputs
    cur_spe <- read_cpout(path, graph_file = "Object_relationships.csv",
                          panel_file = paste0(path, "/panel.csv"))
    
    cur_spe <- read_cpout(path, image_file = NULL, scale_intensities = FALSE, 
                          graph_file = "Object_relationships.csv")
    
    expect_equal(rownames(cur_spe), c("Sm147", "Yb172", "Pr141", "Eu153", "Ag107"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 209))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStack", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    expect_equal(counts(cur_spe), cur_counts)
    expect_equal(counts(cur_spe)[1:10], c(0, 1.41944389703662e-05, 3.54860974259155e-07, 7.09721948518309e-07, 
                                          1.77430487129577e-06, 5.98393015417398e-07, 1.2416655069911e-05, 
                                          1.64558079239784e-06, 8.97589523126097e-07, 1.79517904625219e-06))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`)
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$full, cur_panel$full)
    
    expect_equal(colPairNames(cur_spe), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- vroom::vroom(file.path(path, "Object_relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 3,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_spe$sample_id %in% c("1", "2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "3")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_spe$sample_id %in% c("1", "2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "3")])
    
    expect_silent(cur_graphs <- colPair(cur_spe, "neighborhood"))
    cur_test <- vroom::vroom(file.path(path, "Object_relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 4,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_spe$sample_id %in% c("1", "2", "3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_spe$sample_id == "4")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_spe$sample_id %in% c("1", "2", "3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_spe$sample_id == "4")])
    
    cur_sce <- read_cpout(path, return_as = "sce", image_file = NULL, 
                          scale_intensities = FALSE, graph_file = "Object_relationships.csv")
    
    expect_equal(rownames(cur_sce), c("Sm147", "Yb172", "Pr141", "Eu153", "Ag107"))
    expect_equal(assayNames(cur_sce), "counts")
    expect_equal(dim(cur_sce), c(5, 209))
    expect_equal(names(rowData(cur_sce)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "Pos_X", "Pos_Y", "AreaShape_Area", 
                                            "Neighbors_NumberOfNeighbors_8", "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", 
                                            "AreaShape_MinorAxisLength", "AreaShape_MeanRadius"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStack", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_sce)

    expect_equal(counts(cur_sce), cur_counts)
    expect_equal(counts(cur_sce)[1:10], c(0, 1.41944389703662e-05, 3.54860974259155e-07, 7.09721948518309e-07, 
                                          1.77430487129577e-06, 5.98393015417398e-07, 1.2416655069911e-05, 
                                          1.64558079239784e-06, 8.97589523126097e-07, 1.79517904625219e-06))
    
    expect_equal(cur_sce$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_sce$Pos_X, object_file$Location_Center_X)
    expect_equal(cur_sce$Pos_Y, object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_sce)$Tube.Number, cur_panel$`Tube Number`)
    expect_equal(rowData(cur_sce)$Metal.Tag, cur_panel$`Metal Tag`)
    expect_equal(rowData(cur_sce)$Target, cur_panel$Target)
    expect_equal(rowData(cur_sce)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_sce)$full, cur_panel$full)
    
    expect_equal(colPairNames(cur_sce), "neighborhood")
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- vroom::vroom(file.path(path, "Object_relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 3,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_sce$sample_id %in% c("1", "2")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "3")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_sce$sample_id %in% c("1", "2")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "3")])
    
    expect_silent(cur_graphs <- colPair(cur_sce, "neighborhood"))
    cur_test <- vroom::vroom(file.path(path, "Object_relationships.csv"))
    cur_test <- cur_test[cur_test$`First Image Number` == 4,]
    
    expect_equal(cur_test$`First Object Number` + sum(cur_sce$sample_id %in% c("1", "2", "3")), 
                 from(cur_graphs)[from(cur_graphs) %in% which(cur_sce$sample_id == "4")])
    expect_equal(cur_test$`Second Object Number` + sum(cur_sce$sample_id %in% c("1", "2", "3")), 
                 to(cur_graphs)[to(cur_graphs) %in% which(cur_sce$sample_id == "4")])
    
    cur_spe <- read_cpout(path, scale_intensities = FALSE, graph_file = "Object_relationships.csv")
    
    expect_equal(rownames(cur_spe), c("Sm147", "Yb172", "Pr141", "Eu153", "Ag107"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 209))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius", "ImageNumber", "Metadata_acname", "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStack", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    expect_equal(counts(cur_spe), cur_counts)
    expect_equal(counts(cur_spe)[1:10], c(0, 1.41944389703662e-05, 3.54860974259155e-07, 7.09721948518309e-07, 
                                          1.77430487129577e-06, 5.98393015417398e-07, 1.2416655069911e-05, 
                                          1.64558079239784e-06, 8.97589523126097e-07, 1.79517904625219e-06))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`)
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$full, cur_panel$full)
    
    expect_equal(colPairNames(cur_spe), "neighborhood")
    
    cur_spe <- read_cpout(path, graph_file = NULL)
    
    expect_equal(rownames(cur_spe), c("Sm147", "Yb172", "Pr141", "Eu153", "Ag107"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 209))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius", "ImageNumber", "Metadata_acname", "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStack", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- dplyr::left_join(object_file, image_file, by = "ImageNumber")
    
    expect_equal(counts(cur_spe), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_spe)[1:10], c(0, 0.930232557922948, 0.0232558139480737, 0.0465116278961474, 
                                          0.116279069740369, 0.0392156862653792, 0.813725490006618, 0.107843137229793, 
                                          0.0588235293980688, 0.117647058796138))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_spe$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_spe$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_spe$Metadata_description, object_file$Metadata_description)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`)
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$full, cur_panel$full)
    
    expect_equal(length(colPairNames(cur_spe)), 0)
    
    cur_spe <- read_cpout(path, panel_file = NULL, graph_file = "Object_relationships.csv")
    
    expect_equal(rownames(cur_spe), c("Sm147", "Yb172", "Pr141", "Eu153", "Ag107"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 209))
    expect_equal(length(names(rowData(cur_spe))), 0) 
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius", "ImageNumber", "Metadata_acname", "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStack", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- dplyr::left_join(object_file, image_file, by = "ImageNumber")
    
    expect_equal(counts(cur_spe), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_spe)[1:10], c(0, 0.930232557922948, 0.0232558139480737, 0.0465116278961474, 
                                          0.116279069740369, 0.0392156862653792, 0.813725490006618, 0.107843137229793, 
                                          0.0588235293980688, 0.117647058796138))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_spe$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_spe$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_spe$Metadata_description, object_file$Metadata_description)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    expect_equal(colPairNames(cur_spe), "neighborhood")
    
    cur_spe <- read_cpout(path, intensities = "MedianIntensity_FullStack_", graph_file = "Object_relationships.csv")

    expect_equal(rownames(cur_spe), c("Sm147", "Yb172", "Pr141", "Eu153", "Ag107"))
    expect_equal(assayNames(cur_spe), "counts")
    expect_equal(dim(cur_spe), c(5, 209))
    expect_equal(names(rowData(cur_spe)), 
                 c("Tube.Number", "Metal.Tag", "Target", "ilastik", "full"))
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius", "ImageNumber", "Metadata_acname", "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MedianIntensity_FullStack_", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- merge(object_file, image_file, by = "ImageNumber", order = FALSE)
    
    expect_equal(counts(cur_spe), cur_counts * (2^16 - 1))
    expect_equal(counts(cur_spe)[1:10], c(0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0))
    
    expect_equal(cur_spe$Neighbors_NumberOfNeighbors_8, object_file$Neighbors_NumberOfNeighbors_8)
    expect_equal(cur_spe$Metadata_acname, object_file$Metadata_acname)
    expect_equal(cur_spe$Metadata_acid, object_file$Metadata_acid)
    expect_equal(cur_spe$Metadata_description, object_file$Metadata_description)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,1]), object_file$Location_Center_X)
    expect_equal(as.numeric(spatialCoords(cur_spe)[,2]), object_file$Location_Center_Y)
    
    cur_panel <- vroom::vroom(file.path(path, "panel.csv"))
    mass <- as.numeric(stringr::str_extract(cur_panel$`Metal Tag`, "[0-9]{2,3}"))
    
    expect_equal(rowData(cur_spe)$Tube.Number, cur_panel$`Tube Number`)
    expect_equal(rowData(cur_spe)$Metal.Tag, cur_panel$`Metal Tag`)
    expect_equal(rowData(cur_spe)$Target, cur_panel$Target)
    expect_equal(rowData(cur_spe)$ilastik, cur_panel$ilastik)
    expect_equal(rowData(cur_spe)$full, cur_panel$full)
    
    expect_equal(colPairNames(cur_spe), "neighborhood")
    
    cur_spe <- read_cpout(path, extract_cellmetadata_from = "Location_MaxIntensity_Y_FullStack_c1", 
                          graph_file = "Object_relationships.csv")
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "Location_MaxIntensity_Y_FullStack_c1",
                                            "ImageNumber",
                                            "Metadata_acname", 
                                            "Metadata_acid", "Metadata_description"))
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    
    expect_equal(cur_spe$Location_MaxIntensity_Y_FullStack_c1, object_file$Location_MaxIntensity_Y_FullStack_c1)

    cur_spe <- read_cpout(path, extract_imagemetadata_from = "Metadata_end_timestamp", 
                          graph_file = "Object_relationships.csv")
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius","ImageNumber", "Metadata_end_timestamp"))
    
    expect_equal(spatialCoordsNames(cur_spe), c("Pos_X", "Pos_Y"))
    
    object_file <- vroom::vroom(file.path(path, "cell.csv"))
    cur_counts <- t(object_file[,grepl("MeanIntensity_FullStack", colnames(object_file))])
    rownames(cur_counts) <- rownames(cur_spe)
    
    image_file <- vroom::vroom(file.path(path, "Image.csv"))
    object_file <- merge(object_file, image_file, by = "ImageNumber", order = FALSE)
    
    expect_equal(cur_spe$Metadata_end_timestamp, object_file$Metadata_end_timestamp)
    
    cur_spe <- read_cpout(path, extract_coords_from = NULL, 
                          graph_file = "Object_relationships.csv")
    
    expect_equal(length(spatialCoordsNames(cur_spe)), 0)
    
    cur_sce <- read_cpout(path, extract_coords_from = NULL, return_as = "sce", 
                          graph_file = "Object_relationships.csv")   
    
    expect_equal(names(colData(cur_sce)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius", "ImageNumber", "Metadata_acname", "Metadata_acid", "Metadata_description"
    ))
    
    cur_spe <- read_cpout(path, extract_cellmetadata_from = NULL, 
                          graph_file = "Object_relationships.csv")
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "ImageNumber",
                                            "Metadata_acname", "Metadata_acid", "Metadata_description"))
    
    cur_spe <- read_cpout(path, extract_imagemetadata_from = NULL, 
                          graph_file = "Object_relationships.csv")
    
    expect_equal(names(colData(cur_spe)), c("sample_id", "ObjectNumber", "AreaShape_Area", "Neighbors_NumberOfNeighbors_8", 
                                            "AreaShape_Eccentricity", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength", 
                                            "AreaShape_MeanRadius", "ImageNumber"))
    
    # Fail
    expect_error(read_cpout("test", 
                            graph_file = "Object_relationships.csv"), 
                 regexp = "'path' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(1, 
                            graph_file = "Object_relationships.csv"), 
                 regexp = "'path' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(c("test", "test"), 
                            graph_file = "Object_relationships.csv"), 
                 regexp = "'path' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            object_file = NULL), 
                 regexp = "'object_file' must be specified.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            object_file = "test"), 
                 regexp = "'object_file' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            object_file = 1), 
                 regexp = "'object_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            object_file = c("test", "test")), 
                 regexp = "'object_file' must be a single string.", 
                 fixed = TRUE)

    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            object_feature_file = NULL), 
                 regexp = "'object_feature_file' must be specified.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            object_feature_file = "test"), 
                 regexp = "'object_feature_file' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            object_feature_file = 1), 
                 regexp = "'object_feature_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            object_feature_file = c("test", "test")), 
                 regexp = "'object_feature_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            image_file = NULL), 
                 regexp = "When scaling the summarized object intensities, please supply the 'image_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            image_file = "test"), 
                 regexp = "'image_file' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            image_file = 1), 
                 regexp = "'image_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            image_file = c("test", "test")), 
                 regexp = "'image_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, graph_file = "test"), 
                 regexp = "'graph_file' doesn't exist.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, graph_file = 1), 
                 regexp = "'graph_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, graph_file = c("test", "test")), 
                 regexp = "'graph_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_warning(read_cpout(path, panel_file = "test", 
                              graph_file = "Object_relationships.csv"), 
                 regexp = "'panel_file' does not exist.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            panel_file = 1), 
                 regexp = "'panel_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            panel_file = c("test", "test")), 
                 regexp = "'panel_file' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_metal_from = NULL), 
                 regexp = "'extract_metal_from' must be specified.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_metal_from = "test"), 
                 regexp = "'extract_metal_from' not in panel file.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_metal_from = 1), 
                 regexp = "'extract_metal_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_metal_from = c("test", "test")), 
                 regexp = "'extract_metal_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            intensities =  NULL), 
                 regexp = "'intensities' must be specified.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            intensities =  "test"), 
                 regexp = "No intensity features were read in. Please check the 'intensities' parameter.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            intensities =  "MeanIntensity"), 
                 regexp = "Some of the features set via 'intensities' cannot be uniquely accessed.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            intensities =  c("MeanIntensity", "MedianIntensity")), 
                 regexp = "'intensities' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            intensities =  1), 
                 regexp = "'intensities' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_imgid_from = NULL), 
                 regexp = "'extract_imgid_from' must be specified.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_imgid_from = 1), 
                 regexp = "'extract_imgid_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_imgid_from = c("test", "test")), 
                 regexp = "'extract_imgid_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_imgid_from = "test"), 
                 regexp = "'extract_imgid_from' not in 'object_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_cellid_from = NULL), 
                 regexp = "'extract_cellid_from' must be specified.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_cellid_from = 1), 
                 regexp = "'extract_cellid_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_cellid_from = c("test", "test")), 
                 regexp = "'extract_cellid_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_cellid_from = "test"), 
                 regexp = "'extract_cellid_from' not in 'object_file'.", 
                 fixed = TRUE)

    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_coords_from = 1), 
                 regexp = "'extract_coords_from' not in 'object_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_coords_from = c("test1", "Location_X")), 
                 regexp = "'extract_coords_from' not in 'object_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_cellmetadata_from = 1), 
                 regexp = "'extract_cellmetadata_from' not in 'object_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_cellmetadata_from = c("test1", "Location_X")), 
                 regexp = "'extract_cellmetadata_from' not in 'object_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            scale_intensities = NULL), 
                 regexp = "'scale_intensities' needs to be logical.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            scale_intensities = c(1,2)), 
                 regexp = "'scale_intensities' needs to be of length 1.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            scale_intensities = "test"), 
                 regexp = "'scale_intensities' needs to be logical.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            image_file = NULL), 
                 regexp = "When scaling the summarized object intensities, please supply the 'image_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_imagemetadata_from = "test"), 
                 regexp = "'extract_imagemetadata_from' not in 'image_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_scalingfactor_from = c(1,2)), 
                 regexp = "'extract_scalingfactor_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_scalingfactor_from = "test"), 
                 regexp = "'extract_scalingfactor_from' not in 'image_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_graphimageid_from = NULL), 
                 regexp = "'extract_graphimageid_from' must be specified.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_graphimageid_from = c(1,2)), 
                 regexp = "'extract_graphimageid_from' must be a single string.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_graphimageid_from = "test"), 
                 regexp = "'extract_graphimageid_from' not in 'graph_file'.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_graphcellids_from = NULL), 
                 regexp = "'extract_graphcellids_from' must be specified.", 
                 fixed = TRUE)
    
    expect_error(read_cpout(path, 
                            graph_file = "Object_relationships.csv",
                            extract_graphcellids_from = c("test", "Second Object Number")), 
                 regexp = "'extract_graphcellids_from' not in 'graph_file'.", 
                 fixed = TRUE)
})

