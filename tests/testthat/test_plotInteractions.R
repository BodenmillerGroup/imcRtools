test_that("plotInteractions function works", {
  set.seed(22)
  library(cytomapper)
  data(pancreasSCE)
  
  pancreasSCE <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn",
                                   k = 3)
  
  
  ################################ classic ###################################
  cur_out <- testInteractions(pancreasSCE, 
                              group_by = "ImageNb",
                              label = "CellType",
                              method = "classic",
                              colPairName = "knn_interaction_graph",
                              iter = 100, p_threshold = 0.5,
                              BPPARAM = SerialParam(RNGseed = 123))
  
  ## Plot spatial context - basic tests
  expect_silent(p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb"))
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_size_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_size_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", draw_edges = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "name",
                          node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_color_by = "color_by_sig")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_color_fix = "green")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_width_by = "ct")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_width_fix = 3)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", graph_layout = "linear")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  # return data - tests
  p <- plotInteractions(cur_out,pancreasSCE, label = "CellType",group_by = "ImageNb", return_data = TRUE) 
  expect_type(p, "list")
  expect_equal(names(p), c("edges", "vertices"))
  expect_equal(p$edges[,1], c("celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", "celltype_C", "celltype_C",
                              "celltype_C", "celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", "celltype_C", 
                              "celltype_C", "celltype_C", "celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", 
                              "celltype_C", "celltype_C", "celltype_C"))
  expect_equal(p$edges[,2], c("celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", 
                              "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", 
                              "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", 
                              "celltype_A", "celltype_B", "celltype_C"))
  expect_equal(p$edges[,3], c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(p$vertices[,1], c("celltype_A", "celltype_B", "celltype_C"))
  expect_equal(p$vertices[,2], c(62, 111, 189))
  expect_equal(p$vertices[,3], c(2, 3, 3))
  
  #Errors
  expect_error(plotInteractions(out = cur_out[,1:2], object = pancreasSCE, label = "CellType", group_by = "ImageNb"),
               regexp = "'out' needs to contain columns 'group_by', 'from_label', 'to_label'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = colData(pancreasSCE), label = "CellType", group_by = "ImageNb"),
               regexp = "'object' needs to be a SingleCellExperiment object.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "celltype", group_by = "ImageNb"),
               regexp = "'label' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "imagenb"),
               regexp = "'group_by' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb", 
                                node_color_by = "NAME"),
               regexp = "'node_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_by = "name"),
               regexp = "'node_size_by' has to be 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_fix = factor("black")),
               regexp = "'node_color_fix' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_fix = "3"),
               regexp = "'node_size_fix' has to be numeric.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_repel = "TRUE"),
               regexp = "'node_label_repel' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_by = "NAME"),
               regexp = "'node_label_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_fix = factor("black")),
               regexp = "'node_label_color_fix' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                draw_edges = "TRUE"),
               regexp = "'draw_edges' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                return_data = "TRUE"),
               regexp = "'return_data' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_repel = FALSE, node_label_color_by = "name"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined when node_label_repel == FALSE",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_by = "name", node_color_fix = "blue"),
               regexp = "'node_color_by' and 'node_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_by = "name", node_label_color_fix = "blue"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_by = "name", node_label_color_by = "n_cells"),
               regexp = "'node_label_color_by' and 'node_color_by' have to be identical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_by = "n_group", node_size_fix = 22),
               regexp = "'node_size_by' and 'node_size_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  cur_out[sample(nrow(cur_out), 1), "color_by_sig"] <- NA
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_sig"),
               regexp = "'edge_color_by' contains NA values.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color"),
               regexp = "'edge_color_by' needs to be contained in 'out'.",
               fixed = TRUE)
  
  cur_out$color_by_random <- rownames(cur_out)
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_random"),
               regexp = "'edge_color_by' needs to be unique for all 'from_label'-'to_label' pairs.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_fix = factor("black")),
               regexp = "'edge_color_fix' has to be a character.",
               fixed = TRUE)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_sig", edge_color_fix = "blue"),
               regexp = "'edge_color_by' and 'edge_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = 2),
               regexp = "'edge_width_by' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "weight"),
               regexp = "'edge_width_by' needs to be contained in 'out'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "interaction"),
               regexp = "'edge_width_by' entries need to be numeric.",
               fixed = TRUE)
  
  cur_out$weight <- ifelse(cur_out$from_label == "celltype_A" & cur_out$to_label == "celltype_B", NA, cur_out$ct)
  expect_error(
    plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                     edge_width_by = "weight"),
    regexp = "Missing weights for some 'from_label'-'to_label' pairs",
    fixed = FALSE
  )
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "ct", edge_width_fix = 2),
               regexp = "'edge_width_by' and 'edge_width_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_fix = "3"),
               regexp = "'edge_width_fix' has to be numeric.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                graph_layout = 2),
               regexp = "'graph_layout' has to be a character.",
               fixed = TRUE)

  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                  graph_layout = "unsupported_layout"),
                 regexp = "Layout 'unsupported_layout' not in supported list. Refer to: https://cran.r-project.org/web/packages/ggraph/vignettes/Layouts.html",
                 fixed = TRUE)
  
  
  ################################ conditional ###################################
  cur_out <- testInteractions(pancreasSCE, 
                              group_by = "ImageNb",
                              label = "CellType",
                              method = "conditional",
                              colPairName = "knn_interaction_graph",
                              iter = 100, p_threshold = 0.5,
                              BPPARAM = SerialParam(RNGseed = 123))
  
  ## Plot spatial context - basic tests
  expect_silent(p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb"))
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_size_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_size_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", draw_edges = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "name",
                        node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_color_by = "color_by_sig")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_color_fix = "green")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_width_by = "ct")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_width_fix = 3)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", graph_layout = "linear")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  # return data - tests
  p <- plotInteractions(cur_out,pancreasSCE, label = "CellType",group_by = "ImageNb", return_data = TRUE) 
  expect_type(p, "list")
  expect_equal(names(p), c("edges", "vertices"))
  expect_equal(p$edges[,1], c("celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", "celltype_C", "celltype_C",
                              "celltype_C", "celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", "celltype_C", 
                              "celltype_C", "celltype_C", "celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", 
                              "celltype_C", "celltype_C", "celltype_C"))
  expect_equal(p$edges[,2], c("celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", 
                              "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", 
                              "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", 
                              "celltype_A", "celltype_B", "celltype_C"))
  expect_equal(p$edges[,3], c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(p$vertices[,1], c("celltype_A", "celltype_B", "celltype_C"))
  expect_equal(p$vertices[,2], c(62, 111, 189))
  expect_equal(p$vertices[,3], c(2, 3, 3))
  
  #Errors
  expect_error(plotInteractions(out = cur_out[,1:2], object = pancreasSCE, label = "CellType", group_by = "ImageNb"),
               regexp = "'out' needs to contain columns 'group_by', 'from_label', 'to_label'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = colData(pancreasSCE), label = "CellType", group_by = "ImageNb"),
               regexp = "'object' needs to be a SingleCellExperiment object.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "celltype", group_by = "ImageNb"),
               regexp = "'label' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "imagenb"),
               regexp = "'group_by' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb", 
                                node_color_by = "NAME"),
               regexp = "'node_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_by = "name"),
               regexp = "'node_size_by' has to be 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_fix = factor("black")),
               regexp = "'node_color_fix' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_fix = "3"),
               regexp = "'node_size_fix' has to be numeric.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_repel = "TRUE"),
               regexp = "'node_label_repel' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_by = "NAME"),
               regexp = "'node_label_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_fix = factor("black")),
               regexp = "'node_label_color_fix' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                draw_edges = "TRUE"),
               regexp = "'draw_edges' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                return_data = "TRUE"),
               regexp = "'return_data' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_repel = FALSE, node_label_color_by = "name"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined when node_label_repel == FALSE",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_by = "name", node_color_fix = "blue"),
               regexp = "'node_color_by' and 'node_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_by = "name", node_label_color_fix = "blue"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_by = "name", node_label_color_by = "n_cells"),
               regexp = "'node_label_color_by' and 'node_color_by' have to be identical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_by = "n_group", node_size_fix = 22),
               regexp = "'node_size_by' and 'node_size_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  cur_out[sample(nrow(cur_out), 1), "color_by_sig"] <- NA
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_sig"),
               regexp = "'edge_color_by' contains NA values.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color"),
               regexp = "'edge_color_by' needs to be contained in 'out'.",
               fixed = TRUE)
  
  cur_out$color_by_random <- rownames(cur_out)
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_random"),
               regexp = "'edge_color_by' needs to be unique for all 'from_label'-'to_label' pairs.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_fix = factor("black")),
               regexp = "'edge_color_fix' has to be a character.",
               fixed = TRUE)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_sig", edge_color_fix = "blue"),
               regexp = "'edge_color_by' and 'edge_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = 2),
               regexp = "'edge_width_by' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "weight"),
               regexp = "'edge_width_by' needs to be contained in 'out'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "interaction"),
               regexp = "'edge_width_by' entries need to be numeric.",
               fixed = TRUE)
  
  cur_out$weight <- ifelse(cur_out$from_label == "celltype_A" & cur_out$to_label == "celltype_B", NA, cur_out$ct)
  expect_error(
    plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                     edge_width_by = "weight"),
    regexp = "Missing weights for some 'from_label'-'to_label' pairs",
    fixed = FALSE
  )
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "ct", edge_width_fix = 2),
               regexp = "'edge_width_by' and 'edge_width_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_fix = "3"),
               regexp = "'edge_width_fix' has to be numeric.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                graph_layout = 2),
               regexp = "'graph_layout' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                graph_layout = "unsupported_layout"),
               regexp = "Layout 'unsupported_layout' not in supported list. Refer to: https://cran.r-project.org/web/packages/ggraph/vignettes/Layouts.html",
               fixed = TRUE)
  
  ################################ patch ###################################
  cur_out <- testInteractions(pancreasSCE, 
                              group_by = "ImageNb",
                              label = "CellType",
                              method = "patch",
                              patch_size = 3,
                              colPairName = "knn_interaction_graph",
                              iter = 100, p_threshold = 0.5,
                              BPPARAM = SerialParam(RNGseed = 123))
  
  ## Plot spatial context - basic tests
  expect_silent(p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb"))
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_size_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_size_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", draw_edges = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "name",
                        node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_color_by = "color_by_sig")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_color_fix = "green")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_width_by = "ct")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_width_fix = 3)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", graph_layout = "linear")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  # return data - tests
  p <- plotInteractions(cur_out,pancreasSCE, label = "CellType",group_by = "ImageNb", return_data = TRUE) 
  expect_type(p, "list")
  expect_equal(names(p), c("edges", "vertices"))
  expect_equal(p$edges[,1], c("celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", "celltype_C", "celltype_C",
                              "celltype_C", "celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", "celltype_C", 
                              "celltype_C", "celltype_C", "celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", 
                              "celltype_C", "celltype_C", "celltype_C"))
  expect_equal(p$edges[,2], c("celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", 
                              "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", 
                              "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", 
                              "celltype_A", "celltype_B", "celltype_C"))
  expect_equal(p$edges[,3], c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(p$vertices[,1], c("celltype_A", "celltype_B", "celltype_C"))
  expect_equal(p$vertices[,2], c(62, 111, 189))
  expect_equal(p$vertices[,3], c(2, 3, 3))
  
  #Errors
  expect_error(plotInteractions(out = cur_out[,1:2], object = pancreasSCE, label = "CellType", group_by = "ImageNb"),
               regexp = "'out' needs to contain columns 'group_by', 'from_label', 'to_label'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = colData(pancreasSCE), label = "CellType", group_by = "ImageNb"),
               regexp = "'object' needs to be a SingleCellExperiment object.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "celltype", group_by = "ImageNb"),
               regexp = "'label' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "imagenb"),
               regexp = "'group_by' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb", 
                                node_color_by = "NAME"),
               regexp = "'node_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_by = "name"),
               regexp = "'node_size_by' has to be 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_fix = factor("black")),
               regexp = "'node_color_fix' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_fix = "3"),
               regexp = "'node_size_fix' has to be numeric.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_repel = "TRUE"),
               regexp = "'node_label_repel' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_by = "NAME"),
               regexp = "'node_label_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_fix = factor("black")),
               regexp = "'node_label_color_fix' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                draw_edges = "TRUE"),
               regexp = "'draw_edges' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                return_data = "TRUE"),
               regexp = "'return_data' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_repel = FALSE, node_label_color_by = "name"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined when node_label_repel == FALSE",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_by = "name", node_color_fix = "blue"),
               regexp = "'node_color_by' and 'node_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_by = "name", node_label_color_fix = "blue"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_by = "name", node_label_color_by = "n_cells"),
               regexp = "'node_label_color_by' and 'node_color_by' have to be identical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_by = "n_group", node_size_fix = 22),
               regexp = "'node_size_by' and 'node_size_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  cur_out[sample(nrow(cur_out), 1), "color_by_sig"] <- NA
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_sig"),
               regexp = "'edge_color_by' contains NA values.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color"),
               regexp = "'edge_color_by' needs to be contained in 'out'.",
               fixed = TRUE)
  
  cur_out$color_by_random <- rownames(cur_out)
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_random"),
               regexp = "'edge_color_by' needs to be unique for all 'from_label'-'to_label' pairs.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_fix = factor("black")),
               regexp = "'edge_color_fix' has to be a character.",
               fixed = TRUE)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_sig", edge_color_fix = "blue"),
               regexp = "'edge_color_by' and 'edge_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = 2),
               regexp = "'edge_width_by' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "weight"),
               regexp = "'edge_width_by' needs to be contained in 'out'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "interaction"),
               regexp = "'edge_width_by' entries need to be numeric.",
               fixed = TRUE)
  
  cur_out$weight <- ifelse(cur_out$from_label == "celltype_A" & cur_out$to_label == "celltype_B", NA, cur_out$ct)
  expect_error(
    plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                     edge_width_by = "weight"),
    regexp = "Missing weights for some 'from_label'-'to_label' pairs",
    fixed = FALSE
  )
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "ct", edge_width_fix = 2),
               regexp = "'edge_width_by' and 'edge_width_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_fix = "3"),
               regexp = "'edge_width_fix' has to be numeric.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                graph_layout = 2),
               regexp = "'graph_layout' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                graph_layout = "unsupported_layout"),
               regexp = "Layout 'unsupported_layout' not in supported list. Refer to: https://cran.r-project.org/web/packages/ggraph/vignettes/Layouts.html",
               fixed = TRUE)
  
  ################################ interaction ###################################
  cur_out <- testInteractions(pancreasSCE, 
                              group_by = "ImageNb",
                              label = "CellType",
                              method = "interaction",
                              colPairName = "knn_interaction_graph",
                              iter = 100, p_threshold = 0.5,
                              BPPARAM = SerialParam(RNGseed = 123))
  
  ## Plot spatial context - basic tests
  expect_silent(p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb"))
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_size_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_size_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "name")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "n_group")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_label_color_by = "n_cells")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", draw_edges = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", node_color_by = "name",
                        node_label_repel = FALSE)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_color_by = "color_by_sig")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_color_fix = "green")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_width_by = "ct")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", edge_width_fix = 3)
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  p <- plotInteractions(cur_out, pancreasSCE, "CellType", "ImageNb", graph_layout = "linear")
  expect_s3_class(p, "ggraph")
  expect_silent(print(p))
  expect_equal(p$data$name, sort(unique(pancreasSCE$CellType)))
  expect_equal(p$data$n_cells, as.integer(table(pancreasSCE$CellType)[p$data$name[order(p$data$n_cells)]] %>% unname))
  expect_equal(p$data$n_group, colSums(table(pancreasSCE$ImageNb, pancreasSCE$CellType)[, p$data$name[order(p$data$n_group)]] != 0) %>% unname)
  
  # return data - tests
  p <- plotInteractions(cur_out,pancreasSCE, label = "CellType",group_by = "ImageNb", return_data = TRUE) 
  expect_type(p, "list")
  expect_equal(names(p), c("edges", "vertices"))
  expect_equal(p$edges[,1], c("celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", "celltype_C", "celltype_C",
                              "celltype_C", "celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", "celltype_C", 
                              "celltype_C", "celltype_C", "celltype_A", "celltype_A", "celltype_A", "celltype_B", "celltype_B", "celltype_B", 
                              "celltype_C", "celltype_C", "celltype_C"))
  expect_equal(p$edges[,2], c("celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", 
                              "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", 
                              "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", "celltype_A", "celltype_B", "celltype_C", 
                              "celltype_A", "celltype_B", "celltype_C"))
  expect_equal(p$edges[,3], c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(p$vertices[,1], c("celltype_A", "celltype_B", "celltype_C"))
  expect_equal(p$vertices[,2], c(62, 111, 189))
  expect_equal(p$vertices[,3], c(2, 3, 3))
  
  #Errors
  expect_error(plotInteractions(out = cur_out[,1:2], object = pancreasSCE, label = "CellType", group_by = "ImageNb"),
               regexp = "'out' needs to contain columns 'group_by', 'from_label', 'to_label'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = colData(pancreasSCE), label = "CellType", group_by = "ImageNb"),
               regexp = "'object' needs to be a SingleCellExperiment object.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "celltype", group_by = "ImageNb"),
               regexp = "'label' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "imagenb"),
               regexp = "'group_by' not in 'colData(object)'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb", 
                                node_color_by = "NAME"),
               regexp = "'node_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_by = "name"),
               regexp = "'node_size_by' has to be 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_fix = factor("black")),
               regexp = "'node_color_fix' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_fix = "3"),
               regexp = "'node_size_fix' has to be numeric.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_repel = "TRUE"),
               regexp = "'node_label_repel' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_by = "NAME"),
               regexp = "'node_label_color_by' has to be one off 'name', 'n_cells' or 'n_group'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_fix = factor("black")),
               regexp = "'node_label_color_fix' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                draw_edges = "TRUE"),
               regexp = "'draw_edges' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                return_data = "TRUE"),
               regexp = "'return_data' has to be logical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_repel = FALSE, node_label_color_by = "name"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined when node_label_repel == FALSE",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_by = "name", node_color_fix = "blue"),
               regexp = "'node_color_by' and 'node_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_label_color_by = "name", node_label_color_fix = "blue"),
               regexp = "'node_label_color_by' and 'node_label_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_color_by = "name", node_label_color_by = "n_cells"),
               regexp = "'node_label_color_by' and 'node_color_by' have to be identical.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                node_size_by = "n_group", node_size_fix = 22),
               regexp = "'node_size_by' and 'node_size_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  cur_out[sample(nrow(cur_out), 1), "color_by_sig"] <- NA
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_sig"),
               regexp = "'edge_color_by' contains NA values.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color"),
               regexp = "'edge_color_by' needs to be contained in 'out'.",
               fixed = TRUE)
  
  cur_out$color_by_random <- rownames(cur_out)
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_random"),
               regexp = "'edge_color_by' needs to be unique for all 'from_label'-'to_label' pairs.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_fix = factor("black")),
               regexp = "'edge_color_fix' has to be a character.",
               fixed = TRUE)
  
  cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_color_by = "color_by_sig", edge_color_fix = "blue"),
               regexp = "'edge_color_by' and 'edge_color_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = 2),
               regexp = "'edge_width_by' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "weight"),
               regexp = "'edge_width_by' needs to be contained in 'out'.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "interaction"),
               regexp = "'edge_width_by' entries need to be numeric.",
               fixed = TRUE)
  
  cur_out$weight <- ifelse(cur_out$from_label == "celltype_A" & cur_out$to_label == "celltype_B", NA, cur_out$ct)
  expect_error(
    plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                     edge_width_by = "weight"),
    regexp = "Missing weights for some 'from_label'-'to_label' pairs",
    fixed = FALSE
  )
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_by = "ct", edge_width_fix = 2),
               regexp = "'edge_width_by' and 'edge_width_fix' can not be defined at the same time.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                edge_width_fix = "3"),
               regexp = "'edge_width_fix' has to be numeric.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                graph_layout = 2),
               regexp = "'graph_layout' has to be a character.",
               fixed = TRUE)
  
  expect_error(plotInteractions(out = cur_out, object = pancreasSCE, label = "CellType", group_by = "ImageNb",
                                graph_layout = "unsupported_layout"),
               regexp = "Layout 'unsupported_layout' not in supported list. Refer to: https://cran.r-project.org/web/packages/ggraph/vignettes/Layouts.html",
               fixed = TRUE)
  
}
)
  