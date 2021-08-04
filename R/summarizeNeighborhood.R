#' @rdname summarizeNeighborhood
#' @title Summarizes cell-cell interactions 
#'
#' @description TODO
#'
#' @param object a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#'
#' @section Counting interactions
#' 
#' Explain the three different approaches
#'  
#' @examples 
#' # TODO
#'
#' @author Vito Zanotelli
#' @author Jana Fischer
#' @author adapted by Nils Eling (\email{nils.eling@@uzh.ch})
#'
#' @export
summarizeNeighborhood <- function(object, 
                                 group_by,
                                 label,
                                 method = c("classic", "histocat", "patch"),
                                 patch_size = NULL,
                                 colPairName = NULL){
    
    # Input check
    
    cur_label <- as.factor(colData(object)[[label]])
    cur_table <- .prepare_table(object, group_by, cur_label, colPairName)
    
    # Count interactions
    if (method == "classic") {
        cur_count <- .aggregate_classic(cur_table)
    } else if (method == "histocat") {
        cur_count <- .aggregate_histo(cur_table)
    } else if (method == "patch") {
        cur_count <- .aggregate_classic_patch(cur_table, 
                                              patch_size = patch_size)
    }
    
    return(cur_count)

}
