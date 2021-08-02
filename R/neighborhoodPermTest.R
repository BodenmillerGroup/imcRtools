#' @rdname neighborhoodPermTest
#' @title Performs perturbations to test if two cell-types interact more or less frequently than random.
#'
#' @description This function computes summary statistics (mean, median, and variance) within each factor level.
#'  When multiple entries are provided, summary statistics are calculated for eachcombination of factor entries.
#'
#' @param object a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
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
                                 colPairName = "neighbourhood",
                                 iter = 1000,
                                 min_neighbours = 0){

    # Input check

    # Add neighbouRhood functionality here
    cur_label <- as.factor(colData(object)[[label]])
    cur_table <- .prepare_table(object, img_id, cur_label, colPairName)

}
