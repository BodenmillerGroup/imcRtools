#' @title Visualizes the spatial locations of cells
#'
#' @description A general function to plot spatial locations of cells while 
#' specifying color, shape, size. Cell-cell interactions can be visualized
#' in form of connected points.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object
#' 
#' @return returns a \code{ggplot} object
#' 
#' @examples
#' #TODO
#' 
#' @seealso 
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @export
plotSpatial <- function(object,
                        img_id,
                        coords = c("Pos_X", "Pos_Y"),
                        node_color_by = NULL,
                        node_shape_by = NULL,
                        node_size_by = NULL,
                        edge_color_by = NULL,
                        edge_size_by = NULL,
                        draw_edges = FALSE,
                        arrow = NULL,
                        colPairName = NULL,
                        ncols = NULL,
                        nrows = NULL){
    
    .valid.plotSpatial.input(object, img_id, coords, node_color_by, 
                             node_shape_by, node_size_by, edge_color_by,
                             edge_size_by, draw_edges, arrow, colPairName,
                             ncols, nrows)
    
    nodes <- colData(object)[,c(img_id,node_color_by, node_shape_by, 
                                node_size_by),
                             drop=FALSE]
    nodes[,node_shape_by] <- as.character(nodes[,node_shape_by])
    
    if (draw_edges) {
        edges <- as.data.frame(as(colPair(object, colPairName), "DataFrame"))
        
        if (!is.null(edge_color_by) && 
            edge_color_by %in% colnames(colData(object))) {
            edges[,edge_color_by] <- colData(object)[[edge_color_by]][edges$from]
        }
        
        if (!is.null(edge_size_by) && 
            edge_size_by %in% colnames(colData(object))) {
            edges[,edge_size_by] <- colData(object)[[edge_size_by]][edges$from]
        }
        
        cur_graph <- tbl_graph(nodes = nodes,
                               edges = edges)
    } else {
        cur_graph <- tbl_graph(nodes = nodes)
    }
    
    layout <- create_layout(cur_graph, layout = "manual",
                            x = colData(object)[[coords[1]]],
                            y = colData(object)[[coords[2]]])
    
    # Define column and row number
    if (is.null(ncols) && is.null(nrows)) {
        ni <- length(unique(colData(object)[[img_id]]))
        ncols <- ceiling(sqrt(ni))
        nrows <- ceiling(ncol)
    }
    
    if (draw_edges) {
        if (!is.null(arrow)) {
            p <- ggraph(layout) +
                geom_node_point(aes_string(color = node_color_by,
                                           size = node_size_by,
                                           shape = node_shape_by)) +
                geom_edge_link(aes_string(color = edge_color_by, 
                                          size = edge_size_by),
                               arrow = arrow) +
                facet_nodes(img_id, scale = "free",
                            nrow = nrows, ncol = ncols) 
        } else {
            p <- ggraph(layout) +
                geom_node_point(aes_string(color = node_color_by,
                                           size = node_size_by,
                                           shape = node_shape_by)) +
                geom_edge_link(aes_string(color = edge_color_by, 
                                          size = edge_size_by)) +
                facet_nodes(img_id, scale = "free",
                            nrow = nrows, ncol = ncols) 
        }
    } else {
        p <- ggraph(layout) +
            geom_node_point(aes_string(color = node_color_by,
                                       size = node_size_by,
                                       shape = node_shape_by)) +
            facet_nodes(img_id, scale = "free", 
                        nrow = nrows, ncol = ncols)    
    }
        
    return(p)
}
