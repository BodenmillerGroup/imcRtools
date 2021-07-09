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
                        colPairName = NULL){
    
    #.valid.plotSpatial.input()
    nodes <- colData(object)[,c(img_id,color_by, shape_by, size_by),drop=FALSE]
    nodes[,shape_by] <- as.character(nodes[,shape_by])
    
    if (draw_edges) {
        cur_graph <- tbl_graph(nodes = nodes,
                               edges = as.data.frame(colPair(object, colPairName)))
    } else {
        cur_graph <- tbl_graph(nodes = nodes)
    }
    
    layout <- create_layout(cur_graph, layout = "manual",
                            x = colData(object)[[coords[1]]],
                            y = colData(object)[[coords[2]]])
    
    if (draw_edges) {
        p <- ggraph(layout) +
            geom_node_point(aes_string(color = node_color_by,
                                       size = node_size_by,
                                       shape = node_shape_by)) +
            geom_edge_link(aes()) +
            facet_nodes(img_id, scale = "free") 
    } else {
        p <- ggraph(layout) +
            geom_node_point(aes_string(color = node_color_by,
                                       size = node_size_by,
                                       shape = node_shape_by)) +
            facet_nodes(img_id, scale = "free")    
    }
        
  
    return(p)
}
