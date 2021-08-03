#' @rdname neighborhoodPermTest
#' @title Performs perturbations to test if two cell-types interact more or less frequently than random.
#'
#' @description This function computes summary statistics (mean, median, and variance) within each factor level.
#'  When multiple entries are provided, summary statistics are calculated for eachcombination of factor entries.
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
neighborhoodPermTest <- function(object, 
                                 group_by,
                                 label,
                                 method = c("classic", "histocat", "patch"),
                                 patch_size = NULL,
                                 colPairName = NULL,
                                 iter = 1000,
                                 p_threshold = 0.01,
                                 BBPARAM = SerialParam()){

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
    
    # Permute the labels
    cur_out <- .permute_labels(object, group_by, cur_label, iter,
                               colPairName, method, BBPARAM)
    
    cur_out <- .calc_p_vals(cur_count, cur_out, n_perm = iter, 
                            p_thres = p_threshold)
    
    setorder(cur_out, group_by, from_label, to_label)
    
    return(cur_out)
}
