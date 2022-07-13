#' @title Detect the spatial context of each cell based on its neighborhood
#'
#' @description Function to detect the spatial context (SC) of each cell. 
#' Based on its sorted (high-to-low) cellular neighborhood (CN) fractions in a 
#' k-nearest neighbor graph, the SC of each cell is assigned as the set of CNs 
#' that cumulatively exceed a user-defined fraction threshold. 
#' 
#' The term was coined by \href{https://doi.org/10.1016/j.cels.2021.09.012}{
#' Bhate S. et al., Tissue schematics map the specialization of immune tissue 
#' motifs and their appropriation by tumors, Cell Systems, 2022} and describes 
#' tissue regions in which distinct CNs may be interacting. 
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry single character specifying the \code{colData(object)}
#' entry containing the \code{aggregateNeighbors} DataFrame output. 
#' Defaults to "aggregatedNeighbors". 
#' @param threshold single numeric between 0 and 1 that specifies the fraction 
#' threshold for SC assignment. Defaults to 0.9.
#' @param name single character specifying the name of the output saved in 
#'  \code{colData(object)}. Defaults to "spatial_context".
#'
#' @return returns an object of \code{class(object)} containing a new column 
#' entry to \code{colData(object)[[name]]}
#' 
#' @examples 
#' set.seed(22)
#' library(cytomapper)
#' data(pancreasSCE)
#'
#' ## 1. Cellular neighborhood (CN)
#' sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb",
#'                         type = "knn",
#'                         name = "knn_cn_graph",
#'                         k = 5)
#'
#' sce <- aggregateNeighbors(sce, colPairName = "knn_cn_graph",
#'                          aggregate_by = "metadata",
#'                          count_by = "CellType",
#'                          name = "aggregatedCellTypes")
#'
#' cur_cluster <- kmeans(sce$aggregatedCellTypes, centers = 3)
#' sce$cellular_neighborhood <- factor(cur_cluster$cluster)
#'
#' plotSpatial(sce, img_id = "ImageNb",
#'            colPairName = "knn_cn_graph",
#'            node_color_by = "cellular_neighborhood")
#'
#' ## 2. Spatial context (SC)
#' sce <- buildSpatialGraph(sce, img_id = "ImageNb",
#'                         type = "knn",
#'                         name = "knn_sc_graph",
#'                         k = 15)
#'
#' sce <- aggregateNeighbors(sce, colPairName = "knn_sc_graph",
#'                          aggregate_by = "metadata",
#'                          count_by = "cellular_neighborhood",
#'                          name = "aggregatedNeighborhood")
#'
#' # Detect spatial context
#' sce <- detectSpatialContext(sce, entry = "aggregatedNeighborhood",
#'                            threshold = 0.9)
#'
#' plotSpatial(sce, img_id = "ImageNb",
#'            colPairName = "knn_sc_graph",
#'            node_color_by = "spatial_context")
#'            
#' @seealso 
#' \code{\link[imcRtools]{filterSpatialContext}} for the function to filter 
#' spatial contexts
#'
#' \code{\link[imcRtools]{plotSpatialContext}} for the function to plot 
#' spatial context graphs 
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @references
#' \href{https://doi.org/10.1016/j.cels.2021.09.012}{
#' Bhate S. et al., Tissue schematics map the specialization of immune tissue 
#' motifs and their appropriation by tumors, Cell Systems, 2022}
#' 
#' @importFrom SingleCellExperiment colData
#' @export

detectSpatialContext <- function(object,
                                 entry = "aggregatedNeighbors",
                                 threshold = 0.9,
                                 name = "spatial_context"){
  
  .valid.detectSpatialContext.input(object, entry, threshold, name) 
  
  cur_dat <- colData(object)[,entry]
  
  out_dat <- apply(cur_dat, 1, function(x){
    
    out <- cumsum(sort(x, decreasing = TRUE))
    
    if (sum(out) != 0) {
      cur_out <- names(out[seq_len(sum(out < threshold) + 1)])
      cur_out <- paste(sort(as.numeric(cur_out)), collapse = "_")
      return(cur_out)
    } else {
        return(NA)
    }
  })
  
  colData(object)[[name]] <- out_dat
  return(object)
}