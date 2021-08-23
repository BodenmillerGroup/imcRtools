#' @rdname findBorderCells
#' @title Find cells at the image border
#'
#' @description  
#' Detecting of cells close to the image border for subsequent exlcusion from
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
#' @importFrom data.table setorder
#'
#' @export
findBorderCells <- function(object,
                            img_id,
                            border_dist,
                            coords = c("Pos_X", "Pos_Y"),
                            BPPARAM = SerialParam()){
    
    # Input check
    
    
    return(object)
}
