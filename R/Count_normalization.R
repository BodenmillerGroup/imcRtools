#' @rdname Count_normalization
#' @title Normalisation and scaling of marker expression 
#'
#' @description This function normalize and transform the intensity using a glm-based strategy
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param perform_batch_correction should covariates be regressed out ?
#' @param batch_vector vector describing the batch of each cell
#' @param residual_normalisation method for residual normalisation. Has to be chosen among "Anscombe","Pearson","Working" or "VST"

#' @return Returns a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with a new assay slot called "Normalized_intensity"
#'
#' @examples
#' sce = Count_normalization(sce,residual_normalisation = "Anscombe")
#' @import SingleCellExperiment
#' @import doParallel
#' @import foreach
#' @export

Count_normalization = function(sce,perform_batch_correction=FALSE,
                               batch_vector=NULL,residual_normalisation = "Anscombe") {
  
  #Transforming the data back to count data
  Cell_size = sce$Cell_size
  Transformed_data = t(sce@assays@data[[1]])
  Transformed_data = Transformed_data * 2^metadata(sce)$Bit_mode
  Transformed_data = apply(Transformed_data,MARGIN = 2,FUN = function(x) {x*Cell_size})
  
  Transformed_data = round(Transformed_data)
  
  #Creating parallel backend
  cat(paste("Creating parallel backend using"),as.character(metadata(sce)$N_core),"cores \n")
  registerDoParallel(metadata(sce)$N_core)
  
  #Performing the poisson regression
  if (!perform_batch_correction) {
    cat("Fitting Poisson regressions ...")
    
    List_regression_model =foreach(i=colnames(Transformed_data)) %dopar% {
      Poisson_model = glm(Transformed_data[,i]~log(Cell_size),family = "poisson")
    }
  }
  cat(" done ! \n")
  
  #Extracting and normalizing residuals
  
  
  Residual_matrix =foreach(i=1:length(List_regression_model),.combine = cbind ) %dopar% {
    
    Fitted_values = List_regression_model[[i]]$fitted.values
    Real_values = Transformed_data[,i]
    
    if (!residual_normalisation%in%c("Anscombe","Pearson","Working","VST")){
      cat("No proper method for residual normalization provided. Using the Anscombe normalization method")
    }
    
    if (residual_normalisation=="Anscombe") {
      Normalised_residuals = 1.5 * (Real_values^(2/3)-Fitted_values^(2/3)) / (Fitted_values^(1/6))
    }
    
    if (residual_normalisation=="Pearson") {
      Normalised_residuals = (Real_values-Fitted_values)/sqrt(Fitted_values)
    }
    
    if (residual_normalisation=="Working") {
      Normalised_residuals = (Real_values-Fitted_values)/Fitted_values
    }
    
    if (residual_normalisation=="VST") {
      Normalised_residuals = (sqrt(Real_values)-sqrt(Fitted_values))/2
    }
    
    Normalised_residuals
    
    
  }
  
  #Scaling to 0-1 values
  Residual_matrix = apply(Residual_matrix,MARGIN = 2,FUN = function(x) {
    x = x-min(x)
    x = x/max(x)
  })
  
  Residual_matrix = as.data.frame(Residual_matrix)
  colnames(Residual_matrix) = colnames(Transformed_data)
  
  #Adding the matrix a new assay slot 
  assay(sce, "Normalized_intensity",withDimnames = FALSE) <- t(Residual_matrix)
  
  return(sce)
}

