#' @title Detect the spatial community of each cell
#'
#' @description Function to detect the spatial community of each cell as 
#' proposed by \href{https://www.nature.com/articles/s41586-019-1876-x}{
#' Jackson et al., The single-cell pathology landscape of breast cancer, Nature, 
#' 2020}
#' 
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry single character specifying the \code{colData(object)}
#' entry containing the \code{aggregateNeighbors} DataFrame output. 
#' Defaults to "aggregatedNeighbors". 
#' @param threshold single numeric between 0 and 1 that specifies the fraction 
#' threshold for SC assignment. Defaults to 0.9.
#' @param name single character specifying the name of the output saved in 
#'  \code{colData(object)}. Defaults to "spatial_community".
#'
#' @section Spatial community background:
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
#' plotSpatial(sce, img_id = "ImageNb",
#'            colPairName = "neighborhood",
#'            node_color_by = "spatial_community")
#'            
#' ## Detect spatial community - by celltype
#' sce <- detectCommunity(sce,
#'                       colPairName = "neighborhood,
#'                       group_by = "CellType")
#'
#' plotSpatial(sce, img_id = "ImageNb",
#'            colPairName = "neighborhood",
#'            node_color_by = "spatial_community")
#'            
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @references
#' \href{https://www.nature.com/articles/s41586-019-1876-x}{
#' Jackson et al., The single-cell pathology landscape of breast cancer, 
#' Nature, 2020}
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom stringr str_replace
#' @export


detectCommunity <- function(object,
                            colPairName,
                            size_threshold = 10,
                            group_by = NULL,
                            name = "spatial_community"){
  
#.valid.detectCommunity.input(object, colPairName, size_threshold, group_by, name) 
  
  if (!is.null(group_by)) {
  
  cur_list <- lapply(unique(colData(object)[,group_by]), function(x){
    
    cur_object  <- object[,colData(object)[,group_by] == x]
    gr <- graph_from_data_frame(as.data.frame(colPair(cur_object, colPairName)), 
                                directed = FALSE, 
                                vertices = data.frame(index = seq_len(ncol(cur_object))))
    cl_comm <- cluster_louvain(gr)
    comm <- paste0(x,"_",membership(cl_comm))  
    
    if (size_threshold != 0) {
      comm[membership(cl_comm) %in% which(sizes(cl_comm) < size_threshold)] <- NA
    }
    
    names(comm) <- colnames(cur_object)
    return (comm)
  })
  
  comm <- unlist(cur_list, use.names = TRUE)
  
  } else {
    gr <- graph_from_data_frame(as.data.frame(colPair(object, colPairName)), 
                              directed = FALSE, 
                              vertices = data.frame(index = seq_len(ncol(object))))
    cl_comm <- cluster_louvain(gr)
    comm <- as.character(membership(cl_comm))  
    
    if (size_threshold != 0){
      comm[membership(cl_comm) %in% which(sizes(cl_comm) < size_threshold)] <- NA
      }
    
    names(comm) <- colnames(object)
  }
  
  colData(object)[[name]] <- comm[colnames(object)] 
  
  return(object)
}
