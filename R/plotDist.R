#' @rdname plotDist
#' @title Visualizes distributions of marker intensities
#'
#' @description Ridge plot to visualize the distribution of cell intensities for each marker.
#'  These can be split an colored by cell-level metadata.
#'  By default, the intensity distributions for each marker will be displayed.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param y a character string indicating which factor to plot on the y-axis. Default: y = "rows", which plots a density for each feature.
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

plotDist <- function(x, y = "rows", split_by = NULL, exprs_values = "counts", ...){

  # Check if x is SingleCellExpriment
  .sceCheck(x)
  
  # Check if assay entry exits
  .assayCheck(x, exprs_values)
  
  # The y aesthetic has to be defined
  if(is.null(y)){
    stop("y cannot be empty. Please specify which aesthetic to plot on the y-axis.")
  }

  # Check if selected variable exists
  entries <- c("rows", colnames(colData(x)))
  if(!is.null(y) & !(y %in% entries)){
    stop("The entry for y is not the rownames of colData slots of the object.")
  }

  if(!is.null(split_by) & !(split_by %in% entries)){
    stop("The entry for split_by is not the rownames of colData slots of the object.")
  }

  # Build the data.frame for plotting
  cur_df <- reshape2::melt(assay(x, exprs_values))

  if(y == "rows"){
    cur_df$y <- as.factor(cur_df$Var1)
  }
  else{
    cur_df$y <- as.factor(colData(x)[match(cur_df$Var2, colnames(x)),y])
  }

  if(!is.null(split_by)){
    if(split_by == "rows"){
      cur_df$split <- as.factor(cur_df$Var1)
    }
    else{
      cur_df$split <- as.factor(colData(x)[match(cur_df$Var2, colnames(x)),split_by])
    }

    ggplot(cur_df) + geom_density_ridges(aes_string(x="value", y="y", fill = "y"), ...) + 
      facet_wrap(. ~ split)

  }
  else{
    ggplot(cur_df) + geom_density_ridges(aes_string(x="value", y="y", fill = "y"), ...)

  }

}
