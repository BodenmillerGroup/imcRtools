#' @rdname plotCellCounts
#' @title Visualizes number of cells per factor level
#'
#' @description This function plots the number of cells per chosen cell metadata factor.
#'
#' @param x a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param colour_by character string corresponding to a \code{colData(x)} slot.
#'  The resulting bar plot will be coloured by this factor.
#' @param split_by character string corresponding to a \code{colData(x)} slot.
#'  Defines what will be plotted on the x-axis.
#' @param proportion logical if number of cells per \code{split_by} factor should be scaled to 1.
#'  Default \code{FALSE}.
#'
#' @return Returns a \code{ggplot} barplot object that can be further modified following the \code{ggplot2} syntax.
#'
#' @examples
#' # TODO
#'
#' @author Nils Eling \email{nils.eling@@uzh.ch}
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom SingleCellExperiment colData
#' @importFrom reshape2 melt
#' @export

plotCellCounts <- function(x, colour_by = NULL, split_by = NULL, proportion = FALSE){

  # Check if x is SingleCellExpriment
  .sceCheck(x)

  # Check if selected variable exists
  entries <- colnames(colData(x))
  if(!is.null(colour_by) & !(colour_by %in% entries)){
    stop("The entry for colour_by is not a colData slot of the object.")
  }

  if(!is.null(split_by) & !(split_by %in% entries)){
    stop("The entry for split_by is not a colData slot of the object.")
  }

  # Plot the counts
  if(!is.null(split_by)){
    cur_df <- data.frame(split_by = as.factor(colData(x)[,split_by]),
                         colour_by = as.factor(colData(x)[,colour_by]))
  } else {
    cur_df <- data.frame(split_by = as.factor(rep("All", ncol(x))),
                         colour_by = as.factor(colData(x)[,colour_by]))
  }

  if(proportion){
    df_sum <- melt(table(cur_df)/rowSums(table(cur_df)))

    ggplot(df_sum) + geom_bar(aes(x = split_by, y = value, fill = colour_by), stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            panel.background = element_blank(),
            panel.grid.minor=element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.major.y=element_line(color="grey", size=.3)) +
      ylab("Cell counts") + scale_fill_discrete(name = colour_by)

  } else {
    ggplot(cur_df) + geom_bar(aes(x = split_by, fill = colour_by)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            panel.background = element_blank(),
            panel.grid.minor=element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.major.y=element_line(color="grey", size=.3)) +
      ylab("Cell counts") + scale_fill_discrete(name = colour_by)

  }
}
