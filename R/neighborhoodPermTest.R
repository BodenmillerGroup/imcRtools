#' @rdname neighborhoodPermTest
#' @title Performs perturbations to test if two cell-types interact more or less frequently than random.
#'
#' @description This function computes summary statistics (mean, median, and variance) within each factor level.
#'  When multiple entries are provided, summary statistics are calculated for eachcombination of factor entries.
#'
#' @param object a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param method which aggregation method to use
#'
#' @examples 
#' # TODO
#'
#' @author Vito Zantonelli, adapted by Nils Eling \email{nils.eling@@uzh.ch}
#'
#' @export
neighborhoodPermTest <- function(object, method = c("histocat", "classic")){
    
    stop("This function is under development!")
  # Check if x is SingleCellExpriment
  .sceCheck(object)

  # Add neighbouRhood functionality here

}
