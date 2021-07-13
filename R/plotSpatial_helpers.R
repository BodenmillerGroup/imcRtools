# Function to generate the tidygraph
.generateGraph <- function(object, nodes, colPairName, draw_edges, 
                           edge_color_by, edge_size_by, directed){
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
                               edges = edges,
                               directed = directed)
    } else {
        cur_graph <- tbl_graph(nodes = nodes,
                               directed = directed)
    }
    
    return(cur_graph)
    
}

# Function to generate the base plot
.generatePlot <- function(layout, draw_edges, arrow, node_color_by,
                          node_size_by, node_shape_by, edge_color_by,
                          edge_size_by){
    
    node_color_by <- if(is.null(node_color_by)) NULL else as.name(node_color_by)
    node_size_by <- if(is.null(node_size_by)) NULL else as.name(node_size_by)
    node_shape_by <- if(is.null(node_shape_by)) NULL else as.name(node_shape_by)
    edge_color_by <- if(is.null(edge_color_by)) NULL else as.name(edge_color_by)
    edge_size_by <- if(is.null(edge_size_by)) NULL else as.name(edge_size_by)
    
    if (draw_edges) {
        if (!is.null(arrow)) {
            p <- ggraph(layout) +
                geom_node_point(aes_(colour = node_color_by,
                                     size = node_size_by,
                                     shape = node_shape_by)) +
                geom_edge_link(aes_(edge_colour = edge_color_by, 
                                    size = edge_size_by),
                               arrow = arrow)
        } else {
            p <- ggraph(layout) +
                geom_node_point(aes_(colour = node_color_by,
                                     size = node_size_by,
                                     shape = node_shape_by)) +
                geom_edge_link(aes_(edge_colour = edge_color_by, 
                                    size = edge_size_by)) 
        }
    } else {
        p <- ggraph(layout) +
            geom_node_point(aes_(colour = node_color_by,
                                 size = node_size_by,
                                 shape = node_shape_by))   
    }
    
    return(p)
}
