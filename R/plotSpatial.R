#' @title Visualizes the spatial locations and interactions of cells
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
#' \code{spatialCoords} entries indicating the the cells' x and y locations.
#' @param node_color_by single character indicating the \code{colData(object)}
#' entry or marker name by which the nodes (cell locations) should be colored. 
#' @param node_shape_by single character indicating the \code{colData(object)}
#' entry by which the shape of the nodes are defined. 
#' @param node_size_by single character indicating the \code{colData(object)}
#' entry by which the size of the nodes are defined. 
#' @param node_color_fix single character or numeric specifying the color of all 
#' nodes.
#' @param node_size_fix single numeric specifying the size of all nodes
#' @param node_shape_fix single numeric or character specifying the shape of all 
#' nodes.
#' @param assay_type single character indicating the assay slot from which
#' to extract the expression data when \code{node_color_by} is set to one of
#' \code{rownames(object)}.
#' @param draw_edges should cell-cell interactions be drawn as edges between
#' nodes? 
#' @param directed should cell-cell interactions be handled as a directed graph?
#' @param colPairName single character specifying the 
#' \code{colPair(object)} slot to retrieve the cell-cell pairings.
#' @param edge_color_by single character indicating by which to color the edges.
#' See details for more information.
#' @param edge_width_by single character determining the size of the edges.
#' See details for more information.
#' @param edge_color_fix single character or numeric specifying the color of all 
#' edges.
#' @param edge_width_fix single numeric specifying the size of all edges.
#' @param arrow an \code{\link[grid]{arrow}} object specifying how to draw
#' arrows between cells.
#' @param end_cap a \code{\link[ggraph]{geometry}} object specifying how long
#' the edges are. This only takes effect when drawing arrows. Default:
#' \code{end_cap = circle(0.1, 'cm')}
#' @param nodes_first should the nodes be plotted first and then the edges?
#' @param ncols number of columns of the grid to arrange individual images.
#' @param nrows number of rows of the grid to arrange individual images.
#' @param scales one of \code{"free"}, \code{"fixed"}, \code{"free_x"} or 
#' \code{"free_y"} indicating if x- and y-axis ranges should be fixed across
#' all images. Defaults to "fixed" to match physical units on the x- and y-axis.
#' @param flip_x flip the x-axis?
#' @param flip_y flip the y-axis?
#' @param aspect_ratio relative ratio between the physical units of the x and y 
#' axis (defaults to 1). 
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
#' \code{node_shape_by} and \code{node_size_by} for associating the aesthetics
#' with variables. If node aesthetics should be the same for all nodes,  
#' \code{node_color_fix}, \code{node_shape_fix} and \code{node_size_fix} can 
#' be set.
#' 
#' When \code{draw_edges = TRUE}, cell-cell interactions are visualized in form
#' of edges between nodes. For this, \code{object} needs to contain 
#' column pairings in \code{colPair(object, colPairName)}. Edge color and size
#' can be set by specifying either an entry in 
#' \code{mcols(colPair(object, colPairName))} (edge attributes) or in 
#' \code{colData(object)}. In the latter case, edges are colored by attributes 
#' associated to the "from" node. Variable aesthetics can be set using 
#' \code{edge_color_by} and \code{edge_width_by}. If all edges should have 
#' the same width or color, \code{edge_color_fix} and \code{edge_width_fix}
#' can be set.
#' 
#' Arrows for displaying directed graphs can be drawn by supplying a 
#' \code{\link[grid]{arrow}} object. Arrow attributes can be set within this 
#' class. To cap the edge before it reaches the next node, the \code{end_cap}
#' parameter can be used. 
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
#' # With edges and nodes colored by expression
#' plotSpatial(sce, img_id = "ImageNb",
#'             node_color_by = "PIN",
#'             assay_type = "exprs",
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
#'             edge_color_fix = "green",
#'             arrow = grid::arrow(length = grid::unit(0.1, "inch")),
#'             end_cap = ggraph::circle(0.2, "cm"))
#'   
#' @seealso 
#' \code{\link{buildSpatialGraph}} for constructing interaction graphs
#' 
#' \code{\link{ggraph}} for handling graph aesthetics
#' 
#' @author Nils Eling (\email{nils.eling@@dqbm.uzh.ch})
#' 
#' @import ggraph
#' @importFrom tidygraph tbl_graph  
#' @importFrom ggplot2 aes_ theme element_text element_blank scale_color_manual
#' scale_size_manual scale_shape_manual
#' @export
plotSpatial <- function(object,
                        img_id,
                        coords = c("Pos_X", "Pos_Y"),
                        node_color_by = NULL,
                        node_shape_by = NULL,
                        node_size_by = NULL,
                        node_color_fix = NULL,
                        node_shape_fix = NULL,
                        node_size_fix = NULL,
                        assay_type = NULL,
                        draw_edges = FALSE,
                        directed = TRUE,
                        edge_color_by = NULL,
                        edge_width_by = NULL,
                        edge_color_fix = NULL,
                        edge_width_fix = NULL,
                        arrow = NULL,
                        end_cap = NULL,
                        colPairName = NULL,
                        nodes_first = TRUE,
                        ncols = NULL,
                        nrows = NULL,
                        scales = "fixed",
                        flip_x = FALSE,
                        flip_y = TRUE,
                        aspect_ratio = 1){
    
    .valid.plotSpatial.input(object, img_id, coords, node_color_by, 
                             node_shape_by, node_size_by, edge_color_by,
                             assay_type, edge_width_by, draw_edges, directed, 
                             arrow, end_cap, colPairName, nodes_first,
                             ncols, nrows, scales, flip_x, flip_y, aspect_ratio)
    
    nodes <- .makeNodes(object, node_color_by, img_id, node_shape_by,
                        node_size_by, assay_type)
    
    cur_graph <- .generateGraph(object, nodes, colPairName, draw_edges, 
                                edge_color_by, edge_width_by, directed)
    
    if (is(object, "SpatialExperiment")) {
        layout <- create_layout(cur_graph, layout = "manual",
                                x = spatialCoords(object)[,coords[1]],
                                y = spatialCoords(object)[,coords[2]])
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
    
    p <- .generatePlot(layout, draw_edges, directed, arrow, end_cap, 
                       node_color_by, node_size_by, node_shape_by, 
                       node_color_fix, node_size_fix, node_shape_fix, 
                       edge_color_by, edge_width_by, edge_color_fix, 
                       edge_width_fix, nodes_first)
    
    p <- .postProcessPlot(p, object, img_id, nrows, ncols, node_color_by, 
                          node_color_fix,
                          node_shape_fix, node_size_fix, edge_color_fix, 
                          edge_width_fix, scales, flip_x, flip_y, aspect_ratio)
        
    return(p)
}
