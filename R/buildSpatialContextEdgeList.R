#' @title Create edgelist for spatial context graph
#'
#' @description Function to create a symbolic edge list for spatial context 
#' graph construction. Based on single cell spatial context assignments, this 
#' function operates on cohort-level. 
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry single character specifying the \code{colData(object)} entry 
#' containing the \code{\link[imcRtools]{detectSpatialContext}} output. 
#' Defaults to "spatial_context".
#' @param img_id single character specifying the \code{colData(object)} entry 
#' containing the unique image identifiers. Defaults to "sample_id".
#'
#' @return returns a data frame containing a symbolic edge list in the first
#' two columns ("from" and "to") and, if \code{combined = TRUE}, an additional 
#' column containing the \code{img_id} identifiers.
#' 
#' @examples TO DO
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom stringr str_split
#' @importFrom tibble column_to_rownames
#' @importFrom BiocGenerics table
#' @export

buildSpatialContextEdgeList <- function(object,
                          entry = "spatial_context",
                          img_id = "sample_id")
                          {
  
  .valid.buildEdgeList.input(object, entry, img_id)
  
  data <- colData(object)[,colnames(colData(object)) %in% c(entry,img_id)] %>% table() %>% as.data.frame

  list <- str_split(unique(data$spatial_context), "_")
  list_length <- sapply(list, length)
  edges <- .createEdgeList(list, list_length)
  return(edges)
}