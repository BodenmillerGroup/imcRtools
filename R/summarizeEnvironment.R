#' @title Generates a data frame with the count of cells of a specified type (e.g. clustered celltypes) in the selected neighborhood of each cell.
#'
#' @description Function to generates a data frame with the count of cells of a specified type (e.g. clustered celltypes)
#' in the selected neighborhood of each cell.
#'
#' @param object a \code{SingleCellExperiment} object
#' @param colPairName single character indicating the \code{colPair(object)} entry
#' containing the the neighbor information.
#' @param group single character specifying the \code{colData(object)} entry containing the
#' cellular features that should be summarized in the neighborhood.
#' Supported entries are of type character.
#' @param name single character specifying the name of the data frame.
#'
#' @return returns a \code{SingleCellExperiment}
#' containing the data frame in form of a \code{DataFrame} object in
#' \code{colData(object)}.
#'
#'
#' @author Daniel Schulz (\email{daniel.schulz@@uzh.ch})
#'
#' @importFrom dplyr as_tibble group_by count select
#' @importFrom tidyr pivot_wider
#' @importFrom S4Vectors DataFrame
#' @export

summarizeNeighbors <- function(object,
                               colPairName,
                               group,
                               name = NULL){

  if (! graphName %in% colPairNames(object)) {
    stop("graphName is not a graph in the current object")
  }

  if (! group %in% colnames(colData(object))) {
    stop("group is not a valid colData entry in the current object")
  }

  cur_dat <- as_tibble(colPair(object,graphName))

  cur_out <- cur_dat %>%
    mutate(celltype = colData(object)[[group]][cur_dat$to]) %>%
    group_by(from,.drop=FALSE) %>%
    count(celltype) %>%
    pivot_wider(names_from = celltype,values_from = n, values_fill = 0) %>%
    ungroup() %>%
    select(unique(colData(object)[[group]]))

  # add the dataframe to the object

  name <- ifelse(is.null(name), "summarizedNeighbors", name)

  colData(object)[[name]] <- DataFrame(cur_out)

  return(object)

}
