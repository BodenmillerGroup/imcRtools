#' @title Plot spatial context graph
#'
#' @description Function to plot directed spatial context graphs based on
#' symbolic edge-lists and vertex metadata, which operates on cohort-level.
#' The user can specify node, node_label and edge aesthetics.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object.
#' @param entry single character specifying the \code{colData(object)} entry 
#' containing the \code{detectSpatialContext} output. Defaults to 
#' "spatial_context".
#' @param group_by a single character indicating the \code{colData(object)}
#' entry by which SCs are grouped. This is usually the image or patient ID. 
#' Defaults to "sample_id".
#' @param node_color_by single character from
#' \code{c("name","n_cells","n_group")} by which the nodes should be
#' colored.
#' @param node_size_by single character from \code{c("n_cells","n_group")} 
#' by which the size of the nodes are defined.
#' @param node_color_fix single character specifying the color of all nodes.
#' @param node_size_fix single numeric specifying the size of all nodes.
#' @param node_label_repel should nodes be labelled? Defaults to TRUE.
#' @param node_label_color_by single character from
#' \code{c("name","n_cells","n_group")} by which the node labels should be
#' colored.
#' @param node_label_color_fix single character specifying the color of all node
#' labels.
#' @param draw_edges should edges be drawn between nodes? Defaults to TRUE.
#' @param edge_color_fix single character specifying the color of all edges.
#' @param return_data should the edge list and vertex metadata for graph
#' construction be returned as a \code{list} of two \code{data.frames}?
#' 
#' @return returns a \code{ggplot} object or a \code{list} of two
#' \code{data.frames}.
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
#' # Plot spatial context - default
#' plotSpatialContext(sce, group_by = "ImageNb")
#' 
#' # Plot spatial context - adjust aesthetics
#' plotSpatialContext(sce, group_by = "ImageNb",
#'                    node_color_by = "name",
#'                    node_size_by = "n_cells",
#'                    node_label_color_by = "name")
#'                    
#' plotSpatialContext(sce, group_by = "ImageNb",
#'                    node_color_by = "n_cells",
#'                    node_size_by = "n_group")
#'                    
#' # Plot spatial context - return data
#' plotSpatialContext(sce, group_by = "ImageNb",
#'                   return_data = TRUE)          
#'                   
#' @seealso 
#' \code{\link[imcRtools]{detectSpatialContext}} for the function to detect
#' spatial contexts
#'
#' \code{\link[imcRtools]{filterSpatialContext}} for the function to filter 
#' spatial contexts
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @references
#' \href{https://doi.org/10.1016/j.cels.2021.09.012}{
#' Bhate S. et al., Tissue schematics map the specialization of immune tissue 
#' motifs and their appropriation by tumors, Cell Systems, 2022}
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom stringr str_split
#' @importFrom dplyr count filter group_by_at group_by_all pull select summarise
#' @importFrom igraph graph_from_data_frame
#' @importFrom S4Vectors unfactor
#' @export

plotSpatialContext <- function(object,
                               entry = "spatial_context",
                               group_by = "sample_id",
                               node_color_by = NULL,
                               node_size_by = NULL,
                               node_color_fix = NULL,
                               node_size_fix = NULL,
                               node_label_repel = TRUE,
                               node_label_color_by = NULL,
                               node_label_color_fix = NULL,
                               draw_edges = TRUE,
                               edge_color_fix = NULL,
                               return_data = FALSE){
  
 .valid.plotSpatialContext.input(object, entry, group_by, node_color_by, 
                                  node_size_by, node_color_fix, node_size_fix, 
                                  node_label_repel, node_label_color_by, 
                                  node_label_color_fix, draw_edges, 
                                  edge_color_fix, return_data)
  
 cur_dat <- colData(object) %>% 
   as.data.frame %>% 
   select(all_of(entry), all_of(group_by)) %>% 
   group_by_all() %>% 
   count() %>% 
   as.data.frame
  
  list <- str_split(unique(cur_dat[,entry]), "_")
  list_length <- vapply(list, length, numeric(1))
  edges <- .createEdgeList(list, list_length)
  
  anno <- cur_dat %>% 
    group_by_at(entry) %>% 
    summarise(n_cells = sum(n), n_group = n()) %>% 
    as.data.frame
  
  anno$length <- vapply(str_split(anno[,entry],"_"),length,numeric(1))
  
  graph <- graph_from_data_frame(d = edges, directed = TRUE,vertices = anno)
    
  plot <- .generateSpatialContextPlot(graph = graph, node_color_by, 
                                      node_size_by, node_color_fix, 
                                      node_size_fix, node_label_repel, 
                                      node_label_color_by, node_label_color_fix,
                                      draw_edges, edge_color_fix) 
  
  if (return_data == TRUE) {
      return_data <- list(edges = edges, vertices = anno)
      return(return_data)
  } else { 
      return(plot)
  }  
}