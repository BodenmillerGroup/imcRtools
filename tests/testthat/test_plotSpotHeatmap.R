test_that("plotSpotHeatmap function works.", {
    path <- system.file("extdata/spillover", package = "imcRtools")
    
    # Read in .txt
    expect_silent(cur_sce <- readSCEfromTXT(path, verbose = FALSE))
    
    # Defaults work
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce))
    expect_equal(class(cur_out), "pheatmap")
    
    cur_colours <- matrix(c("#35B779FF", "#424186FF", "#470D60FF", "#440154FF",
                            "#26818EFF", "#FDE725FF", "#20A486FF", "#31688EFF",
                            "#481E70FF", "#3D4D8AFF", "#92D741FF", "#34608DFF",
                            "#3F4788FF", "#355E8DFF", "#25858EFF", "#DDE318FF"),
                          ncol = 4, byrow = TRUE)
    colnames(cur_colours) <- c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")
    rownames(cur_colours) <- c("Dy161", "Dy162", "Dy163", "Dy164")
    expect_equal(cur_out$gtable$grobs[[1]]$children[[1]]$gp$fill,
                 cur_colours)
    
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, log = FALSE))
    cur_colours <- matrix(c("#481E70FF", "#440154FF", "#440154FF", "#440154FF",
                            "#450558FF", "#FDE725FF", "#471063FF", "#440154FF",
                            "#440154FF", "#440154FF", "#355E8DFF", "#440154FF",
                            "#440154FF", "#440154FF", "#450558FF", "#39BA76FF"),
                          ncol = 4, byrow = TRUE)
    colnames(cur_colours) <- c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")
    rownames(cur_colours) <- c("Dy161", "Dy162", "Dy163", "Dy164")
    expect_equal(cur_out$gtable$grobs[[1]]$children[[1]]$gp$fill,
                 cur_colours)
    
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, order_metals = FALSE, cluster_cols = TRUE, cluster_rows = TRUE))
    cur_colours <- matrix(c("#35B779FF", "#440154FF", "#424186FF", "#470D60FF",
                            "#26818EFF", "#31688EFF", "#FDE725FF", "#20A486FF",
                            "#481E70FF", "#34608DFF", "#3D4D8AFF", "#92D741FF",
                            "#3F4788FF", "#DDE318FF", "#355E8DFF", "#25858EFF"),
                          ncol = 4, byrow = TRUE)
    colnames(cur_colours) <- c("Dy161Di", "Dy164Di", "Dy162Di", "Dy163Di")
    rownames(cur_colours) <- c("Dy161", "Dy162", "Dy163", "Dy164")
    expect_equal(cur_out$gtable$grobs[[3]]$children[[1]]$gp$fill,
                 cur_colours)
    
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, threshold = 200))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, threshold = 0))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, log = FALSE, threshold = 200))
    cur_colours <- matrix(c("#FDE725FF", "#440154FF", "#440154FF", "#440154FF",
                            "#FDE725FF", "#FDE725FF", "#FDE725FF", "#440154FF",
                            "#440154FF", "#440154FF", "#FDE725FF", "#440154FF",
                            "#440154FF", "#440154FF", "#FDE725FF", "#FDE725FF"),
                          ncol = 4, byrow = TRUE)
    colnames(cur_colours) <- c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")
    rownames(cur_colours) <- c("Dy161", "Dy162", "Dy163", "Dy164")
    expect_equal(cur_out$gtable$grobs[[1]]$children[[1]]$gp$fill,
                 cur_colours)
    
    expect_silent(cur_out_2 <- plotSpotHeatmap(cur_sce, log = TRUE, threshold = log10(200)))
    expect_equal(cur_out_2$gtable$grobs[[1]]$children[[1]]$gp$fill,
                 cur_colours)
    
    cur_sce$test_spot <- cur_sce$sample_id
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, spot_id = "test_spot"))
    cur_colours <- matrix(c("#35B779FF", "#424186FF", "#470D60FF", "#440154FF",
                            "#26818EFF", "#FDE725FF", "#20A486FF", "#31688EFF",
                            "#481E70FF", "#3D4D8AFF", "#92D741FF", "#34608DFF",
                            "#3F4788FF", "#355E8DFF", "#25858EFF", "#DDE318FF"),
                          ncol = 4, byrow = TRUE)
    colnames(cur_colours) <- c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")
    rownames(cur_colours) <- c("Dy161", "Dy162", "Dy163", "Dy164")
    expect_equal(cur_out$gtable$grobs[[1]]$children[[1]]$gp$fill,
                 cur_colours)
    
    rowData(cur_sce)$test_channel <-  rowData(cur_sce)$channel_name
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, channel_id = "test_channel"))
    expect_equal(cur_out$gtable$grobs[[1]]$children[[1]]$gp$fill,
                 cur_colours)
    
    assay(cur_sce, "exprs") <- log10(counts(cur_sce))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, assay_type = "exprs", log = FALSE))
    
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, statistic = "mean"))
    
    cur_colours <- matrix(c("#39BA76FF", "#404587FF", "#481467FF", "#440154FF",
                            "#277E8EFF", "#FDE725FF", "#20A486FF", "#32658EFF",
                            "#482576FF", "#39558CFF", "#ADDC30FF", "#31688EFF",
                            "#3D4D8AFF", "#34608DFF", "#238A8DFF", "#EBE51AFF"),
                          ncol = 4, byrow = TRUE)
    colnames(cur_colours) <- c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di")
    rownames(cur_colours) <- c("Dy161", "Dy162", "Dy163", "Dy164")
    expect_equal(cur_out$gtable$grobs[[1]]$children[[1]]$gp$fill,
                 cur_colours)
    
    # Passing parameters to pheatmap
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, order_metals = FALSE, 
                                             cluster_rows = TRUE))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, order_metals = FALSE, 
                                             cluster_cols = TRUE))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, color = inferno(50)))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, color = c("red", "green", "blue", "yellow"),
                                             breaks = c(0, 1, 2, 3, 4)))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, color = c("red", "green", "blue", "yellow"),
                                             breaks = c(0, 1, 2, 3, 4), legend_breaks = c(0, 2)))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, cellwidth = 5))
    expect_silent(cur_out <- plotSpotHeatmap(cur_sce, cellwidth = 5, 
                                             annotation_col = data.frame(row.names = c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di"),
                                                                         labels = c("Dy161Di", "Dy162Di", "Dy163Di", "Dy164Di"))))
    
    # Error
    expect_error(plotSpotHeatmap("test"), 
                 regexp = "'object' needs to be a SingleCellExperiment object.", 
                 fixed = TRUE)
    expect_error(plotSpotHeatmap(cur_sce, spot_id = "test"), 
                 regexp = "'spot_id' not in 'colData(object)'.", 
                 fixed = TRUE)
    expect_error(plotSpotHeatmap(cur_sce, channel_id = "test"), 
                 regexp = "'channel_id' not in 'rowData(object)'.", 
                 fixed = TRUE)
    expect_error(plotSpotHeatmap(cur_sce, assay_type = "test"), 
                 regexp = "'assay_type' not in 'assayNames(object)'.", 
                 fixed = TRUE)
    expect_error(plotSpotHeatmap(cur_sce, statistic = "test"), 
                 regexp = "'arg' should be one of “median”, “mean”, “sum”", 
                 fixed = TRUE)
    expect_error(plotSpotHeatmap(cur_sce, log = 1), 
                 regexp = "'log' needs to be logical.", 
                 fixed = TRUE)
    expect_error(plotSpotHeatmap(cur_sce, threshold = "test"), 
                 regexp = "'threshold' needs to be a single numeric.", 
                 fixed = TRUE)
    expect_error(plotSpotHeatmap(cur_sce, order_metals = 1), 
                 regexp = "'order_metals' needs to be logical.", 
                 fixed = TRUE)
    
    
})
