test_that("detectCommunity function works", {
  set.seed(22)
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
  
  sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb",
                         type = "expansion",
                         name = "neighborhood",
                         threshold = 20)
  
  ## Detect spatial community - tests
  # basics
  expect_silent(cur_sce <- detectCommunity(sce, 
                      colPairName = "neighborhood"))
  expect_equal(names(colData(cur_sce)), 
               c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb", 
                 "CellNb", "MaskName", "Pattern", "spatial_community"))
  expect_s4_class(cur_sce , class = "SingleCellExperiment")
  expect_type(cur_sce$spatial_community, "character")
  expect_equal(cur_sce$spatial_community[1:10], c(E34_824 = "1", E34_835 = "2", E34_839 = "1", E34_844 = "1", 
                                                  E34_847 = "2", E34_853 = "1", E34_859 = "1", E34_864 = "1", E34_865 = "1", 
                                                  E34_872 = "1"))
  
  # check that changing size_threshold works
  set.seed(22)
  expect_silent(cur_sce_2 <- detectCommunity(sce, 
                                           colPairName = "neighborhood",
                                           size_threshold = 40))
  expect_false(identical(cur_sce$spatial_community, cur_sce_2$spatial_community))
  expect_equal(cur_sce_2$spatial_community[1:10], c(E34_824 = "1", E34_835 = NA, E34_839 = "1", E34_844 = "1", 
                                                  E34_847 = NA, E34_853 = "1", E34_859 = "1", E34_864 = "1", E34_865 = "1", 
                                                  E34_872 = "1"))
  
  # check that changing name works 
  expect_silent(cur_sce_3 <- detectCommunity(sce, 
                                             colPairName = "neighborhood",
                                             name = "spatial_comm"))
  expect_equal(names(colData(cur_sce_3)), 
               c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb", 
                 "CellNb", "MaskName", "Pattern", "spatial_comm"))
  
  # check that group_by works 
  set.seed(22)
  expect_silent(cur_sce_4 <- detectCommunity(sce, 
                                             colPairName = "neighborhood",
                                             group_by = "CellType",
                                             BPPARAM = SerialParam(RNGseed = 123)
                                             ))
  expect_false(identical(cur_sce_4$spatial_community, cur_sce$spatial_community))
  expect_equal(cur_sce_4$spatial_community[1:10],c(E34_824 = "celltype_C_1", E34_835 = "celltype_C_2", E34_839 = "celltype_C_1", 
                                                   E34_844 = "celltype_C_1", E34_847 = "celltype_C_2", E34_853 = "celltype_C_2", 
                                                   E34_859 = "celltype_C_1", E34_864 = "celltype_C_1", E34_865 = "celltype_C_2", 
                                                   E34_872 = "celltype_C_1"))
  
  
  # manual coding spatial communities gives the same result
  set.seed(22)
  gr <- graph_from_data_frame(as.data.frame(colPair(sce, "neighborhood")), 
                              directed = FALSE, 
                              vertices = data.frame(index = seq_len(ncol(sce))))
  
  clus_comm <- cluster_louvain(gr)
  cl_comm <- as.character(membership(clus_comm))
  cl_comm[membership(clus_comm) %in% which(sizes(clus_comm) < 10)] <- NA
  cur_sce$spatial_community_manual <- cl_comm
  
  expect_equal(unname(cur_sce$spatial_community), cur_sce$spatial_community_manual)
  
  # Different cluster_fun 
  expect_silent(cur_sce <- detectCommunity(cur_sce, 
                                           colPairName = "neighborhood", 
                                           name = "spatial_comm_wt", 
                                           cluster_fun = "walktrap"))
  expect_equal(names(colData(cur_sce)), 
               c("ImageName", "Pos_X", "Pos_Y", "Area", "CellType", "ImageNb", 
                 "CellNb", "MaskName", "Pattern", "spatial_community", 
                 "spatial_community_manual","spatial_comm_wt"))
  
  expect_false(identical(cur_sce$spatial_community, cur_sce$spatial_comm_wt))
  
  # Use different seed for SerialParam()
  set.seed(22)
  expect_silent(cur_sce_5 <- detectCommunity(sce, 
                                             colPairName = "neighborhood",
                                             group_by = "CellType", 
                                             BPPARAM = SerialParam(RNGseed = 22)
                                             ))
  
  expect_false(identical(cur_sce_5$spatial_community, cur_sce_4$spatial_community))
  
  # Errors 
  expect_error(detectCommunity(colData(sce)),
               regexp = "'object' needs to be a SingleCellExperiment object.",
               fixed = TRUE)
  
  expect_error(detectCommunity(sce, colPairName = 10),
               regexp = "'colPairName' must be a single string.",
               fixed = TRUE)
  
  expect_error(detectCommunity(sce, colPairName = "neighbors"),
               regexp = "'colPairName' not in 'colPairNames(object)'.",
               fixed = TRUE)
  
  expect_error(detectCommunity(sce, colPairName = "neighborhood", size_threshold = "10"),
               regexp = "'size_threshold' needs to be a positive single numeric.",
               fixed = TRUE)
  
  expect_error(detectCommunity(sce, colPairName = "neighborhood", group_by = "celltype"),
               regexp = "'group_by' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(detectCommunity(sce, colPairName = "neighborhood", name = c("cool","name")),
               regexp = "'name' has to be a single character.",
               fixed = TRUE)
  
  expect_error(detectCommunity(sce, colPairName = "neighborhood", cluster_fun = c("cool","function")),
               regexp = "'cluster_fun' has to be a single character.",
               fixed = TRUE)
  
  cur_sce_5 <- sce
  colPair(cur_sce_5, "neighborhood") <- colPair(cur_sce_5, "neighborhood")[0,0]
  expect_error(detectCommunity(cur_sce_5, colPairName = "neighborhood"),
               regexp = "No interactions found.",
               fixed = TRUE)
  
  colnames(sce) <- NULL
  expect_error(detectCommunity(sce, colPairName = "neighborhood"),
               regexp = "'colnames' of 'object' need to be specified and unique (e.g. as cell IDs).",
               fixed = TRUE)
}
)