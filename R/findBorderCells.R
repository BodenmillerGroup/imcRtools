#' @rdname findBorderCells
#' @title Find cells at the image border
#'
#' @description  
#' Detecting of cells close to the image border for subsequent exclusion from
#' downstream analyses.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object.
#' @param img_id single character indicating the \code{colData(object)} entry
#' containing the unique image identifiers. 
#' @param border_dist single numeric defining the distance to the image border.
#' The image border here is defined as the minimum and maximum among the 
#' cells' x and y location.
#' @param coords character vector of length 2 specifying the names of the
#' \code{colData} (for a \code{SingleCellExperiment} object) or the
#' \code{spatialCoords} entries of the cells' x and y locations.
#'  
#' @examples 
#' library(cytomapper)
#' data("pancreasSCE")
#' 
#' sce <- findBorderCells(pancreasSCE, img_id = "ImageNb", 
#'                        border_dist = 10)
#' 
#' plotSpatial(sce, img_id = "ImageNb", node_color_by = "border_cells")
#'
#' @author Nils Eling (\email{nils.eling@@uzh.ch})
#' 
#' @importFrom data.table as.data.table setnames
#' @importFrom SpatialExperiment spatialCoords
#'
#' @export
findBorderCells <- function(object,
                            img_id,
                            border_dist,
                            coords = c("Pos_X", "Pos_Y")){
    
    # Input check
    .valid.findBorderCells.input(object, img_id, border_dist, coords)
    
    if (is(object, "SpatialExperiment")) {
        cur_df <- as.data.table(cbind.data.frame(colData(object)[,img_id], 
                                      spatialCoords(object)[,coords]))
    } else {
        cur_df <- as.data.table(colData(object)[,c(img_id, coords)])
    }
    
    Pos_X <- Pos_Y <- border_cells <- NULL
    
    setnames(cur_df, old = names(cur_df), c("img_id", "Pos_X", "Pos_Y")) 
    cur_df[,border_cells := Pos_X <= min(Pos_X) + border_dist | 
                         Pos_X >= max(Pos_X) - border_dist |
                         Pos_Y <= min(Pos_Y) + border_dist |
                         Pos_Y >= max(Pos_Y) - border_dist, by = img_id]
    colData(object)$border_cells <- cur_df$border_cells
    
    return(object)
}
