# Function to generate the tidygraph
.generateGraph <- function(object, nodes, colPairName, draw_edges, 
                           edge_color_by, edge_width_by, directed){
    if (draw_edges) {
        edges <- as.data.frame(as(colPair(object, colPairName), "DataFrame"))
        
        if (!is.null(edge_color_by) && 
            edge_color_by %in% colnames(colData(object))) {
            edges[,edge_color_by] <- colData(object)[[edge_color_by]][edges$from]
        }
        
        if (!is.null(edge_width_by) && 
            edge_width_by %in% colnames(colData(object))) {
            edges[,edge_width_by] <- colData(object)[[edge_width_by]][edges$from]
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
.generatePlot <- function(layout, draw_edges, directed, arrow, end_cap, node_color_by,
                          node_size_by, node_shape_by, node_color_fix, node_size_fix,
                          node_shape_fix, edge_color_by, edge_width_by,
                          edge_color_fix, edge_width_fix){
    
    node_color_by <- if(is.null(node_color_by)) NULL else as.name(node_color_by)
    node_size_by <- if(is.null(node_size_by)) NULL else as.name(node_size_by)
    node_shape_by <- if(is.null(node_shape_by)) NULL else as.name(node_shape_by)
    edge_color_by <- if(is.null(edge_color_by)) NULL else as.name(edge_color_by)
    edge_width_by <- if(is.null(edge_width_by)) NULL else as.name(edge_width_by)
    
    node_color_by <- if(!is.null(node_color_fix)) as.character(node_color_fix) else node_color_by
    node_size_by <- if(!is.null(node_size_fix)) as.character(node_size_fix) else node_size_by
    node_shape_by <- if(!is.null(node_shape_fix)) as.character(node_shape_fix) else node_shape_by
    edge_color_by <- if(!is.null(edge_color_fix)) as.character(edge_color_fix) else edge_color_by
    edge_width_by <- if(!is.null(edge_width_fix)) as.character(edge_width_fix) else edge_width_by
    
    if (draw_edges) {
        if (!is.null(arrow)) { 
            
            if (is.null(end_cap)) {
                end_cap <- circle(0.1, 'cm')
            }
            
            if (directed) {
                cur_geom_edge <- geom_edge_fan(aes_(edge_colour = edge_color_by, 
                                                    edge_width = edge_width_by),
                                               end_cap = end_cap,
                                               arrow = arrow)
            } else {
                cur_geom_edge <- geom_edge_link(aes_(edge_colour = edge_color_by, 
                                                     edge_width = edge_width_by),
                                               end_cap = end_cap,
                                               arrow = arrow)
            }
        } else {
            if (directed) {
                cur_geom_edge <- geom_edge_fan0(aes_(edge_colour = edge_color_by, 
                                                     edge_width = edge_width_by))
            } else {
                cur_geom_edge <- geom_edge_link0(aes_(edge_colour = edge_color_by, 
                                                     edge_width = edge_width_by))
            }
        }
            
        p <- ggraph(layout) + 
             geom_node_point(aes_(colour = node_color_by,
                                  size = node_size_by,
                                  shape = node_shape_by)) +
            cur_geom_edge
    } else {
        p <- ggraph(layout) +
            geom_node_point(aes_(colour = node_color_by,
                                 size = node_size_by,
                                 shape = node_shape_by))   
    }
    
    return(p)
}

# Post process the plots
.postProcessPlot <- function(p, object, img_id, nrows, ncols, node_color_fix,
                             node_shape_fix, node_size_fix, edge_color_fix, 
                             edge_width_fix, scales){
    
    if (!is.null(node_color_fix)) {
        names(node_color_fix) <- as.character(node_color_fix)
        p <- p + scale_color_manual(values = node_color_fix,
                                    guide = "none") 
    }
    if (!is.null(node_shape_fix)) {
        names(node_shape_fix) <- as.character(node_shape_fix)
        p <- p + scale_shape_manual(values = node_shape_fix,
                                    guide = "none") 
    }
    if (!is.null(node_size_fix)) {
        names(node_size_fix) <- as.character(node_size_fix)
        p <- p + scale_size_manual(values = node_size_fix,
                                    guide = "none") 
    }
    if (!is.null(edge_color_fix)) {
        names(edge_color_fix) <- as.character(edge_color_fix)
        p <- p + scale_edge_color_manual(values = edge_color_fix,
                                    guide = "none") 
    }
    if (!is.null(edge_width_fix)) {
        names(edge_width_fix) <- as.character(edge_width_fix)
        p <- p + scale_edge_width_manual(values = edge_width_fix,
                                    guide = "none") 
    }
    
    if (length(unique(colData(object))[[img_id]]) > 1) {
        p <- p + facet_nodes(img_id, scales = scales,
                             nrow = nrows, ncol = ncols) + 
            theme(axis.text = element_text(),
                  panel.background = element_blank()) 
    } else {
        p <- p + theme(axis.text = element_text(),
                  panel.background = element_blank()) +
            ggtitle(unique(colData(object))[[img_id]])
    }
    
}
