test_that("plotInteractionsCozi returns a ggplot object with valid input", {
  df <- data.frame(
    group_by = rep("sample1", 4),
    from_label = rep(c("A", "B"), each = 2),
    to_label = rep(c("X", "Y"), times = 2),
    zscore = c(-2, 0, 1.5, 3),
    cond_ratio = c(0.1, 0.2, 0.3, 0.4),
    sig = c(TRUE, FALSE, TRUE, FALSE)
  )

  p <- plotInteractionsCozi(df)
  expect_s3_class(p, "ggplot")
  expect_true(any(sapply(p$layers, function(l) inherits(l$geom, "GeomPoint"))))
})

test_that("plotInteractionsCozi filters significant interactions if filter_sig=TRUE", {
  df <- data.frame(
    group_by = rep("sample1", 4),
    from_label = rep(c("A", "B"), each = 2),
    to_label = rep(c("X", "Y"), times = 2),
    zscore = c(-2, 0, 1.5, 3),
    cond_ratio = c(0.1, 0.2, 0.3, 0.4),
    sig = c(TRUE, FALSE, TRUE, FALSE)
  )
  
  p <- plotInteractionsCozi(df, filter_sig = TRUE)
  
  plotted_data <- p$data   # the pre-scaled data frame
  expect_equal(plotted_data$cond_ratio, c(0.1, 0.3))
  expect_true(all(plotted_data$sig))  # should all be TRUE now
})

test_that("plotInteractionsCozi respects zscore and dot size limits", {
  df <- data.frame(
    group_by = rep("sample1", 2),
    from_label = c("A", "B"),
    to_label = c("X", "Y"),
    zscore = c(-10, 10),
    cond_ratio = c(0.05, 0.5),
    sig = c(TRUE, TRUE)
  )

  p <- plotInteractionsCozi(
    df,
    zscore_lim = c(-5, 5),
    dot_size_lim = c(0, 1)
  )
  sc <- p$scales$get_scales("size")
  cc <- p$scales$get_scales("colour")

  expect_equal(sc$limits, c(0, 1))
  expect_equal(cc$limits, c(-5, 5))
})

test_that("plotInteractionsCozi errors on missing required columns", {
  bad_df <- data.frame(
    from_label = "A",
    to_label = "B",
    zscore = 1,
    cond_ratio = 0.5
  )
  expect_error(plotInteractionsCozi(bad_df), "missing one or more required columns")
})

