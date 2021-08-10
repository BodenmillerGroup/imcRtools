test_that("summarizeNeighborhood function works", {
    library(cytomapper)
    data(pancreasSCE)
    
    pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                     k = 3)

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

    expect_equal(cur_out$ct, test$ct)    
})
