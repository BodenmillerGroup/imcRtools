#' @rdname Plot_gene_expression_spatial
#' @title Visualization of spatial gene expression
#'
#' @description Plot the spatial distribution of a given gene normalised expression
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param Image_number Number of name of the image/ROI to be plotted
#' @param Gene Name of the gene to be plotted
#'
#' @return Return a plot 
#' @examples 
#' Plot_gene_expression_spatial(sce,Gene = "CD45")
#'
#'
#' @import SingleCellExperiment 


Plot_gene_expression_spatial = function(sce,Image_number = 1,Gene=NULL) {
  if (is.null(Gene) | !Gene%in%rownames(sce@assays@data@listData$Raw_intensity) ) {
    stop("Please select a correct gene to plot ! \n")
  }
  
  if (!"Normalized_intensity"%in%names(assays(sce))) {
    stop("The data have not been normalised. Please proceed to normalisation first")
  }
  
  
  if (Dimension == "2D") {
    Temp_size_data = sqrt(sce$Cell_size[sce$ImageNumber==Image_number])
  }
  
  if (Dimension == "3D") {
    Temp_size_data = (sce$Cell_size[sce$ImageNumber==Image_number])^1/3
  }
  

  Temp_location_data = data.frame(X=sce$Location_Center_X[sce$ImageNumber==Image_number],
                                  Y=sce$Location_Center_Y[sce$ImageNumber==Image_number])
  Temp_expression_data = as.numeric(sce@assays@data@listData$`Normalized_intensity`[Gene,])
  Temp_size_data = sqrt(sce$Cell_size[sce$ImageNumber==Image_number])
  par(las=1,bty="l")
  plot(Temp_location_data,cex=Temp_size_data/10,pch=21,bg=.color_convertion(Temp_expression_data))
  
}
