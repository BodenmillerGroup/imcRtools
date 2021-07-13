#' @title Visualizes the spatial locations of cells
#'
#' @description A general function to plot spatial locations of cells while 
#' specifying color, shape, size. Cell-cell interactions can be visualized
#' in form of edges between points.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object.
#' @param img_id single character indicating the \code{colData(object)} entry
#' containing the unique image identifiers.
#' @param coords character vector of length 2 specifying the names of the
#' \code{colData} (for a \code{SingleCellExperiment} object) or the
#' \code{spatialCoords} entries of the cells' x and y locations.
#' @param node_color_by single character indicating the \code{colData(object)}
#' entry by which the nodes (cell locations) should be colored. 
#' @param node_shape_by single character indicating the \code{colData(object)}
#' entry by which the shape of the nodes are defined. 
#' @param node_size_by single character indicating the \code{colData(object)}
#' entry by which the size of the nodes are defined. 
#' @param draw_edges should cell-cell interactions be drawn as edges between
#' nodes? 
#' @param directed should cell-cell interactions be handled a directed graph?
#' @param colPairName single character specifying from which 
#' \code{colPair(object)} slot to retrieve the cell-cell pairings.
#' @param edge_color_by single character indicating by which to color the edges.
#' See details for more information.
#' @param edge_size_by single character determining the size of the edges.
#' See details for more information.
#' @param arrow a \code{\link[grid]{arrow}} object specifying how to draw arrows
#' between cells.
#' @param ncols number of columns of the grid to arrange individual images.
#' @param nrows number of rows of the grid to arrange individual images.
#' 
#' @return returns a \code{ggplot} object.
#' 
#' @section Visualizing cell locations and cell-cell interactions: 
#' By default, the cells' locations are visualized in form of points (here also
#' referred to as "nodes") on a 2-dimensional plane. The cells' coordinates are
#' extracted either from \code{colData(object)} slot (for a
#' \code{SingleCellExperiment} input object) or from the
#' \code{spatialCoords(object)} slot (for a \code{SpatialExperiment} input
#' object). Node aesthetics are controlled by setting \code{node_color_by}, 
#' \code{node_shape_by} and \code{node_size_by}. 
#' 
#' When \code{draw_edges = TRUE}, cell-cell interactions are visualized in form
#' of edges between nodes. For this, \code{object} needs to contain 
#' column pairings in \code{colPair(object, colPairName)}. Edge color and size
#' can be set by specifying either an entry in 
#' \code{mcols(colPair(object, colPairName))} (edge attributes) or in 
#' \code{colData(object)}. In the latter case, edges are colored by attributes 
#' associated to the "from" node. 
#' 
#' Arrows for displaying directed graphs can be drawn by supplying a 
#' \code{\link[grid]{arrow}} object. Arrow attributes can be set within this 
#' class.
#' 
#' @examples
#' library(cytomapper)
#' data(pancreasSCE)
#' 
#' sce <- buildSpatialGraph(pancreasSCE, img_id = "ImageNb", 
#'                          type = "knn", k = 3, directed = FALSE)
#'
#' # Only nodes
#' plotSpatial(sce, img_id = "ImageNb",
#'             node_color_by = "CellType",
#'             node_shape_by = "ImageNb",
#'             node_size_by = "Area")
#'   
#' # With edges
#' plotSpatial(sce, img_id = "ImageNb",
#'             node_color_by = "CellType",
#'             node_shape_by = "ImageNb",
#'             node_size_by = "Area",
#'             draw_edges = TRUE,
#'             colPairName = "knn_interaction_graph",
#'             edge_color_by = "Pattern")
#'             
#' # With arrows
#' plotSpatial(sce, img_id = "ImageNb",
#'             node_color_by = "CellType",
#'             node_shape_by = "ImageNb",
#'             node_size_by = "Area",
#'             draw_edges = TRUE,
#'             colPairName = "knn_interaction_graph",
#'             edge_color_by = "Pattern",
#'             arrow = grid::arrow(length = unit(0.1, "inch")))
#'   
#' @seealso 
#' \code{\link{buildSpatialGraph}} for constructing interaction graphs
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @importFrom tidygraph tbl_graph 
#' @importFrom ggraph create_layout ggraph geom_node_point geom_edge_link 
#' facet_nodes  
#' @importFrom ggplot2 aes_ theme element_text element_blank 
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
                        directed = TRUE,
                        arrow = NULL,
                        colPairName = NULL,
                        ncols = NULL,
                        nrows = NULL){
    
    .valid.plotSpatial.input(object, img_id, coords, node_color_by, 
                             node_shape_by, node_size_by, edge_color_by,
                             edge_size_by, draw_edges, directed, arrow, colPairName,
                             ncols, nrows)
    
    nodes <- colData(object)[,c(img_id,node_color_by, node_shape_by, 
                                node_size_by),
                             drop=FALSE]
    
    if (!is.null(node_shape_by)) {
        nodes[,node_shape_by] <- as.character(nodes[,node_shape_by])
    }
    
    cur_graph <- .generateGraph(object, colPairName, draw_edges, 
                                edge_color_by, edge_size_by, directed)
    
    if (is(object, "SpatialExperiment")) {
        layout <- create_layout(cur_graph, layout = "manual",
                                x = spatialCoords(object)[[coords[1]]],
                                y = spatialCoords(object)[[coords[2]]])
    } else {
        layout <- create_layout(cur_graph, layout = "manual",
                                x = colData(object)[[coords[1]]],
                                y = colData(object)[[coords[2]]])
    }
    
    # Define column and row number
    if (is.null(ncols) && is.null(nrows)) {
        ni <- length(unique(colData(object)[[img_id]]))
        ncols <- ceiling(sqrt(ni))
        nrows <- ceiling(ncols)
    }
    
    p <- .generatePlot(layout, draw_edges, arrow, node_color_by,
                       node_size_by, node_shape_by, edge_color_by,
                       edge_size_by)
    
    p <- p + facet_nodes(img_id, scales = "free",
                    nrow = nrows, ncol = ncols) + 
        theme(axis.text = element_text(),
              panel.background = element_blank())
        
    return(p)
}
