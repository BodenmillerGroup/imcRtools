#' @rdname findBorderCells
#' @title Find cells at the image border
#'
#' @description  
#' Detecting of cells close to the image border for subsequent exclusion from
#' downstream analyses.
#'
#' @param object a \code{SingleCellExperiment} or \code{SpatialExperiment}
#' object.
#' @param img_id 
#' @param border_dist 
#' @param BPPARAM 
#'  
#' @examples 
#' # TODO
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
                            coords = c("Position_X", "Pos_Y")){
    
    # Input check
    .valid.findBorderCells.input(object, img_id, border_dist, coords)
    
    if (is(cur_obj, "SpatialExperiment")) {
        cur_df <- as.data.table(cbind(colData(object)[,img_id], 
                                      spatialCoords(object)[,coords]))
    } else {
        cur_df <- as.data.table(colData(object)[,c(img_id, coords)])
    }
    
    setnames(cur_df, old = names(cur_df), c("img_id", "Pos_X", "Pos_Y")) 
    cur_df <- cur_df[,border_cells := Pos_X <= min(Pos_X) + border_dist | 
                         Pos_X >= max(Pos_X) - border_dist |
                         Pos_Y <= min(Pos_Y) + border_dist |
                         Pos_Y >= max(Pos_Y) - border_dist, by = img_id]
    colData(object)$border_cells <- cur_df$border_cells
    
    return(object)
}
