#' @title Detect the spatial context of each cell based on its neighborhood
#'
#' @description Function to detect the spatial context (SC) of each cell. 
#' Based on its sorted (high-to-low) cellular neighborhood (CN) fractions in a 
#' k-nearest neighbor graph, the SC of each cell is assigned as the set of CNs 
#' that cumulatively exceed a user-defined fraction threshold. 
#' 
#' The term was coined by Bhates et al. (Cell Systems, 2022) 
#' <https://doi.org/10.1016/j.cels.2021.09.012> and describes tissue regions 
#' in which distinct CNs may be interacting. 
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry single character specifying the \code{colData(object)}
#' entry containing the \code{aggregateNeighbors} DataFrame output. If NULL, 
#' defaults to "aggregatedNeighbors". 
#' @param threshold single numeric between 0 and 1 that specifies the fraction 
#' threshold for SC assignment. Defaults to 0.9.
#' @param name single character specifying the name of the output saved in 
#'  \code{colData(object)}.
#'
#' @return returns an object of \code{class(object)} containing a new column 
#' entry to \code{colData(object)[[name]]}
#' 
#' @examples TO DO
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @importFrom SingleCellExperiment colData
#' @export

detectSpatialContext <- function(object,
                                 entry = NULL,
                                 threshold = 0.9,
                                 name = NULL){
  
  entry <- ifelse(is.null(entry), "aggregatedNeighbors", entry) #default
  name <- ifelse(is.null(name), "spatial_context", name) #default
  
  .valid.detectSpatialContext.input(object, entry, threshold, name) #validity check
  
  cur_dat <- colData(object)[,entry]
  
  out_dat <- apply(cur_dat, 1, function(x){
    
    out <- cumsum(sort(x, decreasing = TRUE))
    
    if(sum(out) != 0){
    return(paste(sort(as.numeric(names(out[seq_len(sum(out < threshold) + 1)]))), collapse = "_"))
    }else{
    return(NA)
    }
  })
  
  colData(object)[[name]] <- out_dat
  return(object)
}