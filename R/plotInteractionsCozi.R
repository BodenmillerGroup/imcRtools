#' @title Plot interaction conditional ratio and z-score
#'
#' @description Function to plot a dot plot visualizing the conditional ratio and
#' z-score of cell-cell interactions calculated with testInteractinos() with method= "conditional"
#' @param out a data frame, usually the output from \code{testInteractions},
#' representing an edge list with columns \code{"group_by", "from_label",
#' "to_label", "zscore", "cond_ratio"}.
#' @param img_id a single character indicating the column in `out` that
#' represents the image or sample ID.
#' @param zscore single character indicating the name of the column in `out`
#' that contains the z-score values. Defaults to "zscore".
#' @param cond_ratio single character indicating the name of the column in
#' `out` that contains the conditional ratio values. Defaults to "cond_ratio".
#' @param filter_sig logical indicating if results should be filtered based on significance testing. 
#' If TRUE, only interactions with sig == TRUE will be plotted. Defaults to FALSE.
#' @param zscore_lim a numeric vector of length 2 specifying the limits for the
#' z-score color scale.
#' @param dot_size_lim a numeric vector of length 2 specifying the limits for the
#' dot size scale.
#'
#' @return returns a \code{ggplot} object.
#'
#' @examples 
#' set.seed(22)
#' library(cytomapper)
#' library(BiocParallel)
#' library(data.table)
#' data(pancreasSCE)
#'
#' ## 1. countInteractions or testInteractions with setting method = "conditional"
#' sce  <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn", k = 3)
#' 
#' count_out <- countInteractions(sce,
#'                                group_by = "ImageNb",
#'                                label = "CellType",
#'                                method = "conditional", 
#'                                colPairName = "knn_interaction_graph")
#' 
#' test_out <- testInteractions(sce, 
#'                              group_by = "ImageNb",
#'                              label = "CellType", 
#'                              method = "conditional", 
#'                              colPairName = "knn_interaction_graph", 
#'                              iter = 100, 
#'                              p_threshold = 0.5, 
#'                              BPPARAM = SerialParam(RNGseed = 123))
#' 
#' ## 2. Plot interactions as dotplot
#' 
#' # default                
#' plotInteractionsCozi(test_out)
#' 
#' # adjust zscore and dot size limits
#' plotInteractionsCozi(test_out,
#'                     zscore_lim = c(-50, 50),
#'                    dot_size_lim = c(0, 0.5))  
#' 
#' # filter for significant interactions only
#' plotInteractionsCozi(test_out,
#'                    filter_sig = TRUE)   
#'                   
#' @seealso 
#' \code{\link{countInteractions}} for counting (but not testing) cell-cell
#' interactions per grouping level.
#' \code{\link{testInteractions}} for testing cell-cell 
#' interactions per grouping level.
#'
#' @author Chiara Schiller (\email{chiara.schiller@uni-heidelberg.de})
#'
#' @importFrom RColorBrewer brewer.pal
#' @export

plotInteractionsCozi <- function(out,
                                 img_id = "group_by",
                                 zscore = "zscore",
                                 cond_ratio = "cond_ratio",
                                 sig = "sig",
                                 filter_sig = FALSE,
                                 zscore_lim = NULL,
                                 dot_size_lim = NULL) {

                                  # Call the validity check function first
    .valid.plotInteractionsCozi.input(out,
                                      img_id,
                                      zscore,
                                      cond_ratio,
                                      sig,
                                      filter_sig,
                                      zscore_lim,
                                      dot_size_lim)

    required_cols <- c(img_id, "from_label", "to_label", zscore, cond_ratio)
    if (filter_sig) {
        required_cols <- c(required_cols, "sig")
    }
    if (!all(required_cols %in% colnames(out))) {
        stop("Input data frame is missing one or more required columns. Make sure you ran count and testInteractions() with method = 'conditional'.")
    }
    out = as.data.frame(out)

    plot_data <- out %>%
        select(all_of(img_id),
               from_label,
               to_label,
               sig = all_of(sig),
               zscore = all_of(zscore),
               cond_ratio = all_of(cond_ratio))
    
    if (filter_sig) {
        plot_data <- plot_data[plot_data$sig == TRUE, ]
    }

    plot_data <- na.omit(plot_data)
    
    p <- ggplot(plot_data, aes(x = from_label, y = to_label)) +
        geom_point(aes(size = cond_ratio, color = zscore), shape = 16) +
        
        scale_size_continuous(range = c(1, 10), limits = dot_size_lim) +
        ggplot2::scale_color_gradient2(
            low = RColorBrewer::brewer.pal(11, "RdBu")[10], # Use the first color for 'low'
            mid = "white",
            high = RColorBrewer::brewer.pal(11, "RdBu")[2], # Use the last color for 'high'
            midpoint = 0,
            limits = zscore_lim,
            oob = scales::squish
        ) +
        
        facet_wrap(~get(img_id)) +

        labs(
            title = "COZI Z-score and Conditional Ratio",
            x = "Cell Type (From)",
            y = "Cell Type (To)",
            size = "Conditional Ratio",
            color = "Z-score"
        ) +

        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
        )

    return(p)
}