# Helper functions to check entries and validity of SingleCellExperiment
#' @importFrom methods is
.sceCheck <- function(x){
  if(!is(x, "SingleCellExperiment")){
    stop("x is not a SingleCellExperiment object.")
  }
}

.assayCheck <- function(x, exprs_values){
  .sceCheck(x)
  if(!(exprs_values %in% names(assays(x)))){
    stop(paste("The", exprs_values, "slot does not exists in the SingleCellExperiment object."))
  }
}
