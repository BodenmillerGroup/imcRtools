test_that("summarizeNeighborhood function works", {
    library(cytomapper)
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    ########################### classic ############################

    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "CellType", method = "classic",
                                                   colPairName = "knn_interaction_graph")) 
    cur_out <- as.data.frame(cur_out)
    cur_out <- cur_out[order(cur_out$group_by, cur_out$from_label, cur_out$to_label),]
    
    # Check against manual calculations
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["CellType"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["CellType"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(n))
    
    test <- as.data.frame(test)
    test <- test[order(test$image, test$from_label, test$`as.factor(to_label)`),]

    expect_equal(cur_out$ct[!is.na(cur_out$ct)], test$ct[!is.na(cur_out$ct)])   

    # As factor
    pancreasSCE$CellType <- as.factor(pancreasSCE$CellType)
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "CellType", method = "classic",
                                                   colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
        
    expect_equal(cur_out_2$ct, cur_out$ct)
    
    # As character
    pancreasSCE$CellType <- as.character(pancreasSCE$CellType)
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "CellType", method = "classic",
                                                     colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    expect_equal(cur_out_2, cur_out)
    
    # As numeric
    pancreasSCE$CellType <- as.numeric(as.factor(pancreasSCE$CellType))
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "CellType", method = "classic",
                                                     colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    expect_equal(cur_out_2$ct, cur_out$ct)
    
    # Logical
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "Pattern", method = "classic",
                                                   colPairName = "knn_interaction_graph")) 

    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["Pattern"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["Pattern"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(n))
    
    expect_equal(cur_out_2$ct[!is.na(cur_out_2$ct)], test$ct[!is.na(cur_out_2$ct)]) 
    
    # One image only contains one cell type
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    pancreasSCE$test <- pancreasSCE$CellType
    pancreasSCE$test[pancreasSCE$ImageNb == 3] <- "test" 
    expect_silent(cur_out_3 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "test", method = "classic",
                                                   colPairName = "knn_interaction_graph")) 
    
    cur_out_3 <- as.data.frame(cur_out_3)
    cur_out_3 <- cur_out_3[order(cur_out_3$group_by, cur_out_3$from_label, cur_out_3$to_label),]
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["test"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["test"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(n))
    
    expect_equal(cur_out_3$ct[!is.na(cur_out_3$ct)], test$ct[!is.na(cur_out_3$ct)]) 
    expect_equal(cur_out_3$ct[!is.na(cur_out_3$ct) & cur_out_3$group_by != 3], cur_out$ct[!is.na(cur_out$ct) & cur_out$group_by != 3]) 
    
    data(pancreasSCE)
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    pancreasSCE$ImageName <- "test"
    
    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageName",
                                                   label = "CellType", method = "classic",
                                                   colPairName = "knn_interaction_graph")) 
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["CellType"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["CellType"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageName"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(n))
    
    expect_equal(cur_out$ct, test$ct)
    
    ########################### histocat ############################
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "CellType", method = "histocat",
                                                   colPairName = "knn_interaction_graph")) 
    
    cur_out <- as.data.frame(cur_out)
    cur_out <- cur_out[order(cur_out$group_by, cur_out$from_label, cur_out$to_label),]
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["CellType"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["CellType"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label)) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(n))
    
    test <- as.data.frame(test)
    test <- test[order(test$image, test$from_label, test$`as.factor(to_label)`),]
    
    expect_equal(cur_out$ct, test$ct)
    
    # As factor
    pancreasSCE$CellType <- as.factor(pancreasSCE$CellType)
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "CellType", method = "histocat",
                                                     colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    expect_equal(cur_out_2$ct, cur_out$ct)
    
    # As character
    pancreasSCE$CellType <- as.character(pancreasSCE$CellType)
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "CellType", method = "histocat",
                                                     colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    expect_equal(cur_out_2, cur_out)
    
    # As numeric
    pancreasSCE$CellType <- as.numeric(as.factor(pancreasSCE$CellType))
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "CellType", method = "histocat",
                                                     colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    expect_equal(cur_out_2$ct, cur_out$ct)
    
    # Logical
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "Pattern", method = "histocat",
                                                   colPairName = "knn_interaction_graph")) 
    
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["Pattern"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["Pattern"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label)) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(n))
    
    expect_equal(cur_out_2$ct, test$ct) 
    
    # One image only contains one cell type
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    pancreasSCE$test <- pancreasSCE$CellType
    pancreasSCE$test[pancreasSCE$ImageNb == 3] <- "test" 
    expect_silent(cur_out_3 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "test", method = "histocat",
                                                     colPairName = "knn_interaction_graph")) 
    
    cur_out_3 <- as.data.frame(cur_out_3)
    cur_out_3 <- cur_out_3[order(cur_out_3$group_by, cur_out_3$from_label, cur_out_3$to_label),]
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["test"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["test"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label)) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(n))
    
    expect_equal(cur_out_3$ct, test$ct) 
    expect_equal(cur_out_3$ct[cur_out_3$group_by != 3], cur_out$ct[cur_out$group_by != 3]) 
    
    data(pancreasSCE)
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    pancreasSCE$ImageName <- "test"
    
    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageName",
                                                   label = "CellType", method = "histocat",
                                                   colPairName = "knn_interaction_graph")) 
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["CellType"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["CellType"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageName"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label)) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(n))
    
    expect_equal(cur_out$ct, test$ct)
    
    ########################### patch ############################
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "CellType", method = "patch",
                                                   patch_size = 3,
                                                   colPairName = "knn_interaction_graph")) 
    
    cur_out <- as.data.frame(cur_out)
    cur_out <- cur_out[order(cur_out$group_by, cur_out$from_label, cur_out$to_label),]
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["CellType"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["CellType"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>% ungroup() %>% 
        mutate(ct = n >= 3) %>%
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(ct))
    
    test <- as.data.frame(test)
    test <- test[order(test$image, test$from_label, test$`as.factor(to_label)`),]
    
    expect_equal(cur_out$ct[!is.na(cur_out$ct)], test$ct[!is.na(cur_out$ct)])
    
    # As factor
    pancreasSCE$CellType <- as.factor(pancreasSCE$CellType)
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "CellType", method = "patch",
                                                     patch_size = 3,
                                                     colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    expect_equal(cur_out_2$ct, cur_out$ct)
    
    # As character
    pancreasSCE$CellType <- as.character(pancreasSCE$CellType)
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "CellType", method = "patch",
                                                     patch_size = 3,
                                                     colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    expect_equal(cur_out_2, cur_out)
    
    # As numeric
    pancreasSCE$CellType <- as.numeric(as.factor(pancreasSCE$CellType))
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "CellType", method = "patch",
                                                     patch_size = 3,
                                                     colPairName = "knn_interaction_graph")) 
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    expect_equal(cur_out_2$ct, cur_out$ct)
    
    # Logical
    expect_silent(cur_out_2 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "Pattern", method = "patch",
                                                     patch_size = 3,
                                                     colPairName = "knn_interaction_graph")) 
    
    cur_out_2 <- as.data.frame(cur_out_2)
    cur_out_2 <- cur_out_2[order(cur_out_2$group_by, cur_out_2$from_label, cur_out_2$to_label),]
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["Pattern"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["Pattern"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>%
        mutate(ct = n >= 3) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(ct))
    
    expect_equal(cur_out_2$ct[!is.na(cur_out_2$ct)], test$ct[!is.na(cur_out_2$ct)]) 
    
    # One image only contains one cell type
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    pancreasSCE$test <- pancreasSCE$CellType
    pancreasSCE$test[pancreasSCE$ImageNb == 3] <- "test" 
    expect_silent(cur_out_3 <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                     label = "test", method = "patch",
                                                     patch_size = 3,
                                                     colPairName = "knn_interaction_graph")) 
    
    cur_out_3 <- as.data.frame(cur_out_3)
    cur_out_3 <- cur_out_3[order(cur_out_3$group_by, cur_out_3$from_label, cur_out_3$to_label),]
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["test"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["test"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>% 
        mutate(ct = n >= 3) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(ct))
    
    expect_equal(cur_out_3$ct[!is.na(cur_out_3$ct)], test$ct[!is.na(cur_out_3$ct)]) 
    expect_equal(cur_out_3$ct[!is.na(cur_out_3$ct) & cur_out_3$group_by != 3], cur_out$ct[!is.na(cur_out$ct) & cur_out$group_by != 3]) 
    
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    
    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "CellType", method = "patch",
                                                   patch_size = 4,
                                                   colPairName = "knn_interaction_graph")) 
    
    expect_true(all(cur_out$ct[!is.na(cur_out$ct)] == 0))
    
    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "CellType", method = "patch",
                                                   patch_size = 0,
                                                   colPairName = "knn_interaction_graph")) 
    
    expect_true(all(cur_out$ct[!is.na(cur_out$ct)] == 1))
    
    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageNb",
                                                   label = "CellType", method = "patch",
                                                   patch_size = 1,
                                                   colPairName = "knn_interaction_graph")) 
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["CellType"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["CellType"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageNb"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>% 
        mutate(ct = n >= 1) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(ct))
    
    expect_equal(cur_out$ct[!is.na(cur_out$ct)], test$ct[!is.na(cur_out$ct)]) 
    
    data(pancreasSCE)
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)
    pancreasSCE$ImageName <- "test"
    
    expect_silent(cur_out <- summarizeNeighborhood(pancreasSCE, group_by = "ImageName",
                                                   label = "CellType", method = "patch",
                                                   patch_size = 3,
                                                   colPairName = "knn_interaction_graph")) 
    
    cur_table <- as.data.frame(colPair(pancreasSCE, "knn_interaction_graph"))
    cur_table$from_label <- colData(pancreasSCE)[["CellType"]][cur_table$from]
    cur_table$to_label <- colData(pancreasSCE)[["CellType"]][cur_table$to]
    cur_table$image <- colData(pancreasSCE)[["ImageName"]][cur_table$from]
    
    test <- cur_table %>% group_by(from, image, from_label) %>%
        count(as.factor(to_label), .drop = FALSE) %>% 
        mutate(ct = n >= 3) %>% ungroup() %>% 
        group_by(image, from_label, `as.factor(to_label)`) %>% summarize(ct = mean(ct))
    
    expect_equal(cur_out$ct, test$ct)
    
    # Fail
    expect_error(summarizeNeighborhood("test"),
                 regexp = "'object' not of type 'SingleCellExperiment'.",
                 fixed = TRUE)
    expect_error(summarizeNeighborhood(pancreasSCE, group_by = "test", label = "CellType", colPairName = "knn_interaction_graph"),
                 regexp = "'group_by' not in colData(object).",
                 fixed = TRUE)
    expect_error(summarizeNeighborhood(pancreasSCE, group_by = 1, label = "CellType", colPairName = "knn_interaction_graph"),
                 regexp = "'group_by' must be a single string.",
                 fixed = TRUE)
    expect_error(summarizeNeighborhood(pancreasSCE, group_by = "ImageNb", label = "test", colPairName = "knn_interaction_graph"),
                 regexp = "'label' not in colData(object).",
                 fixed = TRUE)
    expect_error(summarizeNeighborhood(pancreasSCE, group_by = "ImageNb", label = 1, colPairName = "knn_interaction_graph"),
                 regexp = "'label' must be a single string.",
                 fixed = TRUE)
    expect_error(summarizeNeighborhood(pancreasSCE, group_by = "ImageNb", label = "CellType", colPairName = "test"),
                 regexp = "'colPairName' not in colPairNames(object).",
                 fixed = TRUE)
    expect_error(summarizeNeighborhood(pancreasSCE, group_by = "ImageNb", label = "CellType", colPairName = 1),
                 regexp = "'colPairName' must be a single string.",
                 fixed = TRUE)
    expect_error(summarizeNeighborhood(pancreasSCE, group_by = "ImageNb", label = "CellType", colPairName = "knn_interaction_graph",
                                       method = "patch"),
                 regexp = "When method = 'patch', please specify 'patch_size'.",
                 fixed = TRUE)
    expect_error(summarizeNeighborhood(pancreasSCE, group_by = "ImageNb", label = "CellType", colPairName = "knn_interaction_graph",
                                       method = "patch", patch_size = "test"),
                 regexp = "'patch_size' must be a single numeric.",
                 fixed = TRUE)
    
})
