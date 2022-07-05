#' @title Plot spatial context graph
#'
#' @description Function to plot directed spatial context graphs based on
#' symbolic edge-lists and vertex metadata, which operates on cohort-level.
#' The user can specify node, node_label and edge aesthetics.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' @param entry single character specifying the \code{colData(object)} entry 
#' containing the \code{detectSpatialContext} output. Defaults to 
#' "spatial_context".
#' @param img_id single character specifying the \code{colData(object)} entry 
#' containing the unique image identifiers. Defaults to "sample_id".
#' @param node_color_by single character from
#' \code{c("name","n_cells","n_samples")} by which the nodes should be
#' colored.
#' @param node_size_by single character from \code{c("n_cells","n_samples")} 
#' by which the size of the nodes are defined.
#' @param node_color_fix single character specifying the color of all nodes.
#' @param node_size_fix single numeric specifying the size of all nodes.
#' @param node_label_repel should nodes be labelled? Defaults to TRUE.
#' @param node_label_color_by single character from
#' \code{c("name","n_cells","n_samples")} by which the node labels should be
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
#' @examples TO DO
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom stringr str_split
#' @importFrom dplyr count filter group_by_at pull select summarise
#' @importFrom igraph graph_from_data_frame
#' @importFrom BiocGenerics table
#' @export

plotSpatialContext <- function(object,
                               entry = "spatial_context",
                               img_id = "sample_id",
                               #nodes
                               node_color_by = NULL,
                               node_size_by = NULL,
                               node_color_fix = NULL,
                               node_size_fix = NULL,
                               #node labels
                               node_label_repel = TRUE,
                               node_label_color_by = NULL,
                               node_label_color_fix = NULL,
                               #plot graph - edges
                               draw_edges = TRUE,
                               edge_color_fix = NULL,
                               #return data
                               return_data = FALSE){
  
  .valid.plotSpatialContext.input(object, entry, img_id, node_color_by, 
                                  node_size_by, node_color_fix, node_size_fix, 
                                  node_label_repel, node_label_color_by, 
                                  node_label_color_fix, draw_edges, 
                                  edge_color_fix, return_data)

  data <- colData(object)[,colnames(colData(object)) %in% c(entry,img_id)] %>% 
  table() %>% as.data.frame
  Freq <- as.name("Freq")
  n <- as.name("n")
  
  list <- str_split(unique(data[,entry]), "_")
  list_length <- sapply(list, length)
  edges <- .createEdgeList(list, list_length)
  
  anno <- data.frame(spatial_context = unique(data[,entry]) %>% unfactor(), 
                     length = listLen(str_split(unique(data[,entry]),"_")),
                     n_cells = data %>% group_by_at(entry) %>% 
                       summarise(sum = sum(Freq)) %>% pull(sum), 
                     n_samples = data %>% group_by_at(entry) %>% 
                       filter(Freq != 0) %>% count() %>% pull(n) %>% 
                       as.integer()
                     )
  
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
