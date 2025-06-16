#' @title Plot interaction graph
#'
#' @description Function to plot directed interaction graphs based on symbolic 
#' edge-lists and vertex metadata.
#' The user can specify node, node_label and edge aesthetics.
#'
#' @param out a data frame representing an edge list with columns \code{"group_by",
#' "from_label" and "to_label"}. Additional columns may be included to specify 
#' edge attributes (weight or color). 
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object.
#' @param label single character specifying the \code{colData(object)} entry
#' which stores the cell labels. These can be cell-types labels or other
#' metadata entries.
#' @param group_by a single character indicating the \code{colData(object)}
#' entry by which interactions are grouped. This is usually the image or patient ID. 
#' @param node_color_by single character either
#' \code{NULL, "name","n_cells", "n_group"} by which the nodes should be
#' colored.
#' @param node_size_by single character either \code{NULL, "n_cells","n_group"} 
#' by which the size of the nodes are defined.
#' @param node_color_fix single character specifying the color of all nodes.
#' @param node_size_fix single numeric specifying the size of all nodes.
#' @param node_label_repel should nodes be labelled? Defaults to TRUE.
#' @param node_label_color_by single character either
#' \code{NULL, "name","n_cells","n_group"} by which the node labels should be
#' colored.
#' @param node_label_color_fix single character specifying the color of all node
#' labels.
#' @param edge_color_by single character indicating the column name of \code{"out"} 
#' by which the edges are colored.
#' @param edge_color_fix single character specifying the color of all edges.
#' @param edge_width_by single character indicating the column name of \code{"out"} 
#' by which the width of the edges are scaled.
#' @param edge_width_fix single numeric specifying the width of all edges.
#' @param draw_edges should edges be drawn between nodes? Defaults to TRUE.
#' @param graph_layout single character either
#' \code{"circle", "chord", "linear", "fr", "kk", "drl", "stress", "graphopt", 
#' "lgl", "tree", "sugiyama", "star", "nicely", "manual", "grid", "mds", "sphere", 
#' "randomly", "gem", "dh"} which defines the graph layout.
#' Defaults to \code{"circle"}. For more information, see \link[ggraph]{ggraph}.
#' @param return_data should the edge list and vertex metadata for graph
#' construction be returned as a \code{list} of two \code{data.frames}?
#' 
#' @return returns a \code{ggplot} object or a \code{list} of two
#' \code{data.frames}.
#' 
#' @examples 
#' set.seed(22)
#' library(cytomapper)
#' library(BiocParallel)
#' library(dplyr)
#' library(ggplot2)
#' library(ggraph)
#' data(pancreasSCE)
#'
#' ## 1. countInteractions or testInteractions
#' sce  <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", type = "knn", k = 3)
#' 
#' cur_out <- countInteractions(sce,
#'                              group_by = "ImageNb",
#'                              label = "CellType",
#'                              method = "classic", # choose from c("classic", "histocat", "patch", "interaction")
#'                              colPairName = "knn_interaction_graph")
#' 
#' cur_out <- testInteractions(sce, 
#'                             group_by = "ImageNb",
#'                             label = "CellType", 
#'                             method = "classic", # choose from c("classic", "histocat", "patch", "interaction")
#'                             patch_size = 3,
#'                             colPairName = "knn_interaction_graph", 
#'                             iter = 100, 
#'                             p_threshold = 0.5, 
#'                             BPPARAM = SerialParam(RNGseed = 123))
#' 
#' ## 2. Plot interactions 
#' 
#' # default                
#' plotInteractions(cur_out, sce, "CellType", "ImageNb")
#' 
#' # adjust node aesthetics
#' plotInteractions(cur_out, sce, "CellType", "ImageNb",
#'                  node_color_by = "name",
#'                  node_size_by = "n_cells")
#'                  
#' ## specify custom node and node label colors
#' node_color <- setNames(c("red", "blue", "green"), c("celltype_A", "celltype_B", "celltype_C"))
#' plotInteractions(cur_out, sce, "CellType", "ImageNb", node_color_by = "name", node_label_color_by = "name") +
#' scale_color_manual(values = node_color)
#'                  
#' # adjust edge aesthetics
#' ## edge width based on the summarized count derived from \code{countInteractions} averaged per 'from_label'-'to_label' pair
#' plotInteractions(cur_out, sce, "CellType", "ImageNb",
#'                  edge_width_by = "ct")
#'                  
#' ## edge color based on interaction significance derived from \code{testInteractions}
#' cur_out <- cur_out %>% as.data.frame() %>% group_by(from_label, to_label) %>% mutate(color_by_sig = sum(sigval, na.rm = TRUE) > 0)
#' plotInteractions(cur_out, sce, "CellType", "ImageNb",
#'                  edge_color_by = "color_by_sig")
#'                  
#' ## specify custom edge colors
#' edge_color <- setNames(c("red", "blue"),c(TRUE, FALSE))
#' plotInteractions(cur_out, sce, "CellType", "ImageNb",
#'                  edge_color_by = "color_by_sig") +
#'                  scale_edge_colour_manual(values = edge_color)
#'                    
#' # Plot spatial context - return data
#' plotInteractions(cur_out, sce, "CellType", "ImageNb",
#'                  return_data = TRUE)          
#'                   
#' @seealso 
#' #' \code{\link{countInteractions}} for counting (but not testing) cell-cell
#' interactions per grouping level.
#' #' \code{\link{testInteractions}} for testing cell-cell 
#' interactions per grouping level.
#' 
#' @author Lasse Meyer (\email{lasse.meyer@@uzh.ch})
#' 
#' @references
#' \href{https://doi.org/10.1016/j.cels.2021.09.012}{
#' Bhate S. et al., Tissue schematics map the specialization of immune tissue 
#' motifs and their appropriation by tumors, Cell Systems, 2022}
#' 
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr %>% group_by summarise filter mutate select n count left_join across group_by_at
#' @importFrom tidyselect all_of
#' @importFrom tidyr separate
#' @importFrom igraph graph_from_data_frame
#' @importFrom stats na.omit
#' @export

plotInteractions <- function(out,
                             object,
                             label,
                             group_by,
                             node_color_by = NULL,
                             node_size_by = NULL,
                             node_color_fix = NULL,
                             node_size_fix = NULL,
                             node_label_repel = TRUE,
                             node_label_color_by = NULL, 
                             node_label_color_fix = NULL,
                             edge_color_by = NULL,
                             edge_color_fix = NULL,
                             edge_width_by = NULL,
                             edge_width_fix = NULL,
                             draw_edges = TRUE,
                             return_data = FALSE,
                             graph_layout = "circle") {
  
  .valid.plotInteractions.input(out, object, label, group_by, node_color_by, node_size_by,
                                node_color_fix, node_size_fix, node_label_repel,
                                node_label_color_by, node_label_color_fix,
                                edge_color_by, edge_color_fix, edge_width_by, edge_width_fix,
                                draw_edges, graph_layout, return_data)
  
  # Edge attributes
  if (is.null(edge_width_by)) {
    edges <- out %>% as.data.frame() %>%
      select(from_label, to_label, color = all_of(edge_color_by)) %>%
      mutate(weight = 1)
  } else {
    edges <- out %>% as.data.frame() %>%
      select(from_label, to_label, color = all_of(edge_color_by), weight = all_of(edge_width_by)) %>%
      group_by(from_label, to_label) %>%
      mutate(weight = mean(weight, na.rm = TRUE))
  }
  
  # Node attributes
  cur_dat <- colData(object) %>%
    as.data.frame() %>%
    select(all_of(label), all_of(group_by)) %>%
    group_by(across(everything())) %>%
    count() %>%
    na.omit() %>%
    as.data.frame()
  
  anno <- cur_dat %>%
    group_by_at(label) %>%
    summarise(n_cells = sum(n), n_group = n(), .groups = "drop") %>%
    filter(.data[[label]] %in% c(edges$from_label, edges$to_label)) %>%
    as.data.frame
  
  # Create graph object
  graph <- graph_from_data_frame(d = edges, directed = TRUE, vertices = anno)
  
  # Generate plot
  plot <- .generateInteractionsPlot(graph = graph,
                                    node_color_by = node_color_by,
                                    node_size_by = node_size_by,
                                    node_color_fix = node_color_fix,
                                    node_size_fix = node_size_fix,
                                    node_label_repel = node_label_repel,
                                    node_label_color_by = node_label_color_by,
                                    node_label_color_fix = node_label_color_fix,
                                    edge_color_by = edge_color_by,
                                    edge_color_fix = edge_color_fix,
                                    edge_width_by = edge_width_by,
                                    edge_width_fix = edge_width_fix,
                                    draw_edges = draw_edges,
                                    graph_layout = graph_layout)
  
  # Return data or plot
  if (return_data) {
    return(list(edges = edges, vertices = anno))
  } else {
    return(plot)
  }
}
