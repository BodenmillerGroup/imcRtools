#' @title Detect the spatial community of each cell
#'
#' @description Function to detect the spatial community of each cell as 
#' proposed by \href{https://www.nature.com/articles/s41586-019-1876-x}{
#' Jackson et al., The single-cell pathology landscape of breast cancer, Nature, 
#' 2020}. Each cell is clustered (Louvain community detection) based on its 
#' interactions as defined by a spatial object graph.
#' 
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param colPairName single character indicating the \code{colPair(object)}
#' entry containing the neighbor information. 
#' @param size_threshold single positive numeric that specifies the minimum 
#' number of cells per community. Defaults to 10.
#' @param group_by single character indicating that spatial community 
#' detection will be performed separately for all unique entries to 
#' \code{colData(object)[,group_by]}.
#' @param name single character specifying the name of the output saved in 
#'  \code{colData(object)}. Defaults to "spatial_community".
#'
#' @section Spatial community detection procedure:
#' 1. Create an igraph object from the edge list stored in 
#' \code{colPair(object, colPairName)}
#' 2. Perform community detection using Louvain algorithm
#' 3. Store the community IDs in a vector and replace all communities with 
#' a size smaller than \code{size_threshold} by NA
#' 
#' Optional steps: Specify \code{group_by} to perform spatial community 
#' detection separately for all unique entries to \code{colData(object)[,group_by]} 
#' e.g. for tumor and stromal cells.
#' 
#' @return returns an object of \code{class(object)} containing a new column 
#' entry to \code{colData(object)[[name]]}
#' 
#' @examples 
#' set.seed(22)
#' library(cytomapper)
#' data(pancreasSCE)
#'
#' sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb",
#'                         type = "expansion",
#'                         name = "neighborhood",
#'                         threshold = 20)
#'
#' ## Detect spatial community 
#' sce <- detectCommunity(sce, 
#'                       colPairName = "neighborhood")
#'
#' plotSpatial(sce,
#'             img_id = "ImageNb",
#'             colPairName = "neighborhood",
#'             node_color_by = "spatial_community")
#'            
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @references
#' \href{https://www.nature.com/articles/s41586-019-1876-x}{
#' Jackson et al., The single-cell pathology landscape of breast cancer, 
#' Nature, 2020}
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom igraph membership sizes
#' @export


detectCommunity <- function(object,
                            colPairName,
                            size_threshold = 10,
                            group_by = NULL,
                            name = "spatial_community"){
  
.valid.detectCommunity.input(object, colPairName, size_threshold, group_by, name) 
  
  if (!is.null(group_by)) {
    cur_list <- lapply(unique(colData(object)[,group_by]), function(x){
      cur_object  <- object[,colData(object)[,group_by] == x]
      cl_comm <- .detectCommunity_function(cur_object, colPairName)
      comm <- paste0(x,"_",membership(cl_comm))  
      names(comm) <- colnames(cur_object)
      
      if (size_threshold != 0){
        comm[membership(cl_comm) %in% which(sizes(cl_comm) < size_threshold)] <- NA
      }
      
      return(comm)
  })
    comm <- unlist(cur_list, use.names = TRUE)
  
  } else {
    cur_object <- object
    cl_comm <- .detectCommunity_function(cur_object, colPairName)
    comm <- as.character(membership(cl_comm))  
    names(comm) <- colnames(object)
    
    if (size_threshold != 0){
      comm[membership(cl_comm) %in% which(sizes(cl_comm) < size_threshold)] <- NA
    }
  }
  
  colData(object)[[name]] <- comm[colnames(object)] 
  
  return(object)
}
