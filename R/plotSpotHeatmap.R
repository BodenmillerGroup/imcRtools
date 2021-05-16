#' @title Summarizes and visualizes the pixel intensities per spot and channel
#'
#' @description Helper function for estimating the spillover matrix. This
#' function visualizes the median pixel intensities per spot (rows) and per
#' channel (columns) in form of a heatmap.
#'
#' @param object a \code{SingleCellExperiment} object containing pixel
#' intensities per channel. Individual pixels are stored as columns and
#' channels are stored as rows.
#' @param spot_id character string indicating which \code{colData(object)} entry
#' stores the isotope names of the spotted metal. Entries should be of the 
#' form (mt)(mass) (e.g. Sm152 for Samarium isotope with the atomic mass 152).
#' @param channel_id character string indicating which \code{rowData(object)} 
#' entry contains the isotope names of the acquired channels. 
#' @param assay_type character string indicating which assay to use (default
#' \code{counts}).
#' @param statistic the statistic to use when aggregating channels per spot
#' (default \code{median})
#' @param log should the aggregated pixel intensities be \code{log10(x + 1)} 
#' transformed?
#' @param threshold single numeric indicating a threshold after pixel
#' aggregation. All aggregated values larger than \code{threshold} will be
#' labeled as \code{1}.
#' @param order_metals should the metals be ordered based on spotted mass?
#' @param color see parameter in \code{\link[pheatmap]{pheatmap}}
#' @param breaks see parameter in \code{\link[pheatmap]{pheatmap}}
#' @param legend_breaks see parameter in \code{\link[pheatmap]{pheatmap}}
#' @param cluster_cols see parameter in \code{\link[pheatmap]{pheatmap}}
#' @param cluster_rows see parameter in \code{\link[pheatmap]{pheatmap}}
#' @param ... other arguments passed to \code{pheatmap}.
#' 
#' @section Quality control for spillover estimation:
#' Visualizing the aggregated pixel intensities serves two purposes:
#' 
#' \enumerate{
#' \item Small median pixel intensities (< 200 counts) might hinder the robust
#' estimation of the channel spillover. In that case, consecutive pixels can be
#' summed (see \code{\link{binAcrossPixels}}).
#' \item Each spotted metal (row) should show the highest median pixel intensity
#' in its corresponding channel (column). If this is not the case, either the
#' naming of the .txt files was incorrect or the incorrect metal was spotted.
#' }
#' 
#' By setting the \code{threshold} parameter, the user can easily identify spots
#' where pixel intensities are too low for robust spillover estimation. 
#' 
#' @return a \code{\link[pheatmap]{pheatmap}} object
#'
#' @examples
#' path <- system.file("extdata/spillover", package = "imcRtools")
#' # Read in .txt files
#' sce <- readSCEfromTXT(path)
#' 
#' # Visualizes heatmap
#' plotSpotHeatmap(sce)
#' 
#' # Visualizes thresholding results
#' plotSpotHeatmap(sce, log = FALSE, threshold = 200)
#' 
#' @seealso \code{\link[pheatmap]{pheatmap}} for visual modifications
#' @seealso \code{\link[scuttle]{aggregateAcrossCells}} for the aggregation
#' function
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#'
#' @importFrom pheatmap pheatmap
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom viridis viridis
#' @importFrom stringr str_extract
#' @importFrom SummarizedExperiment assay
#' @export
plotSpotHeatmap <- function(object, 
                            spot_id = "sample_id",
                            channel_id = "channel_name",
                            assay_type = "counts",
                            statistic = c("median", "mean", "sum"),
                            log = TRUE,
                            threshold = NULL,
                            order_metals = TRUE,
                            color = viridis(100),
                            breaks = NA,
                            legend_breaks = NA,
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            ...){
    
    .valid.plotSpotHeatmap.input(object, spot_id, channel_id, assay_type, 
                           log, threshold, order_metals)
    
    statistic <- match.arg(statistic)
    
    cur_out <- aggregateAcrossCells(object, object[[spot_id]], 
                                    statistics = statistic,
                                    use.assay.type = assay_type)
    
    if (log) {
        cur_mat <- log10(assay(cur_out, assay_type) + 1)
    } else {
        cur_mat <- assay(cur_out, assay_type)
    }
    
    if (!is.null(threshold)) {
        cur_mat <- (cur_mat > threshold) * 1
        
        if (is.na(breaks)) { 
            breaks <- c(0, 0.5, 1)
        }
        
        if (is.na(legend_breaks)) { 
            legend_breaks <- c(0, 1)
        }

        color <- c(color[1], color[length(color)])
    }
    
    colnames(cur_mat) <- cur_out[[spot_id]] 
    rownames(cur_mat) <- rowData(cur_out)[[channel_id]] 
    
    # Order rows and cols based on spot metal
    if (order_metals) {
        cur_spots <- colnames(cur_mat)
        cur_mass <- as.numeric(str_extract(cur_spots, "[0-9]{2,3}$"))
        cur_spots <- cur_spots[order(cur_mass)]
        
        cur_channels <- rownames(cur_mat)
        cur_mass <- as.numeric(str_extract(cur_channels, "[0-9]{2,3}"))
        cur_channels <- cur_channels[order(cur_mass)]
        cur_isotope <- str_extract(cur_channels, "[A-Za-z]{1,2}[0-9]{2,3}")
        cur_rownames <- c(cur_channels[cur_isotope %in% cur_spots],
                          cur_channels[!cur_isotope %in% cur_spots])
        
        cur_mat <- cur_mat[cur_rownames,cur_spots]
        
        cluster_cols <- FALSE
        cluster_rows <- FALSE
    }
    
    # Transposed to match the CATALYST visualization
    pheatmap(t(cur_mat), color = color, 
            cluster_cols = cluster_cols,
            cluster_rows = cluster_rows, 
            breaks = breaks, legend_breaks = legend_breaks,
            ...)
    
}
