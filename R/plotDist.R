#' @rdname plotDist
#' @title Visualizes distributions of marker intensities
#'
#' @description Ridge plot to visualize the distribution of cell intensities for each marker.
#'  These can be split an colored by cell-level metadata.
#'  By default, the intensity distributions for each marker will be displayed.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param colour_by a character string indicating which factor to use for separating the counts.
#'  Default: colour_by = "rows", which plots a density for each feature.
#'  Valid arguments are "rows" or any entry in the \code{colData(sce)} slot.
#' @param split_by character string corresponding to a \code{colData(x)} or \code{color_by = "rows"}, calling the rownames for \code{x}.
#'  This feature will be used to facet wrap the plots.
#' @param exprs_values character string indicating from which \code{assays(x)} slot the intensity values can be extracted.
#'  Default: exprs_values = "counts".
#' @param ... further parameters for the \code{\link[ggridges]{ggridges}} function.
#'
#' @return Returns a \code{ggplot} object that can be further modified following the \code{ggplot2} syntax.
#'
#' @examples
#' # TODO
#'
#' @author Nils Eling \email{nils.eling@@uzh.ch}
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment colData
#' @importFrom reshape2 melt
#' @export

plotDist <- function(x, colour_by = "rows", split_by = NULL,
                     exprs_values = "counts", plot_type = c("ridges", "boxplot"), ...){

  # Check if x is SingleCellExpriment
  .sceCheck(x)

  # Check if assay entry exits
  .assayCheck(x, exprs_values)

  # Select plot type
  plot_type <- match.arg(plot_type)

  # The colour_by aesthetic has to be defined
  if(is.null(colour_by)){
    stop("colour_by cannot be empty. Please specify by which colData entry the counts are separated.")
  }

  # Check if selected variable exists
  entries <- c("rows", colnames(colData(x)))
  if(!is.null(colour_by) & !(colour_by %in% entries)){
    stop("The entry for colour_by is not the rownames of colData slots of the object.")
  }

  # Build the data.frame for plotting
  cur_df <- reshape2::melt(assay(x, exprs_values))

  if(colour_by == "rows"){
    cur_df$colour_by <- as.factor(cur_df$Var1)
  }
  else{
    cur_df$colour_by <- as.factor(colData(x)[match(cur_df$Var2, colnames(x)),colour_by])
  }

  if(!is.null(split_by)){

    if(!(split_by %in% entries)){
      stop("The entry for split_by is not the rownames of colData slots of the object.")
    }

    if(split_by == "rows"){
      cur_df$split <- as.factor(cur_df$Var1)
    }
    else{
      cur_df$split <- as.factor(colData(x)[match(cur_df$Var2, colnames(x)),split_by])
    }

    if(plot_type == "ridges"){
      ggplot(cur_df) +
        geom_density_ridges(aes_string(x="value", y="colour_by", fill = "colour_by"), ...) +
        facet_wrap(. ~ split) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
    }
    else{
      ggplot(cur_df) +
        geom_boxplot(aes_string(x="colour_by", y="value", fill = "colour_by"), ...) +
        facet_wrap(. ~ split) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
    }

  }
  else{
    if(plot_type == "ridges"){
      ggplot(cur_df) +
        geom_density_ridges(aes_string(x="value", y="colour_by", fill = "colour_by"), ...) +
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
    }
    else{
      ggplot(cur_df) +
        geom_boxplot(aes_string(x="colour_by", y="value", fill = "colour_by"), ...) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }

}
