#' @rdname Create_SCE
#' @title Generates an SCE data container for multiplexed single-cell imaging data 
#' @description This function takes the list produced by the Load_IMC_data() function and outputs an SCE object
#' @param List_data a list object produced by Load_IMC_data
#' @param dimension the dimension of the data, either 2D or 3D
#' @param Bit_mode bit mode acquisition, usually set to 16
#' @param N_core number of cores used during parallelized computational steps

#' @return returns an SCE object
#'
#' @examples
#' sce = Create_SCE(List_data,dimension = "2D",Bit_mode = 16,N_core = 6)
#'
#' @import SingleCellExperiment
#' @import S4Vectors
#' @export

Create_SCE = function(List_data,dimension = "2D",Bit_mode=16,N_core = 6) {
  
  sce = SingleCellExperiment(assays = list(Raw_intensity = as.matrix(t(List_data$Expression_data))),
                             colData = List_data$Cell_annotation,
                             rowData = List_data$Gene_annotation,
                             metadata = list(dimension = dimension,Bit_mode=Bit_mode,N_core = N_core))
  return(sce)
}

