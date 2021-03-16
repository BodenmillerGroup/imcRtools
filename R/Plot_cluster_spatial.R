#' @rdname Plot_cluster_spatial
#' @title Visualization of cluster spatial distribution
#'
#' @description Plot the spatial distribution of a given gene normalised expression
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Clustering of the cell should have been previously performed.
#' @param Image_number Number of name of the image/ROI to be plotted
#'
#' @return Return a plot 
#' @examples 
#' Plot_cluster_spatial(sce,Image_number = 1)
#'
#'
#' @import SingleCellExperiment 


Plot_cluster_spatial = function(sce,Image_number = 1) {
  if (is.null(colLabels(sce))) {
    stop("Please compute clustering first ! \n")
  }
  
  
  Temp_location_data = data.frame(X=sce$Location_Center_X[sce$ImageNumber==Image_number],
                                  Y=sce$Location_Center_Y[sce$ImageNumber==Image_number])
  Temp_cluster_data = colLabels(sce)[sce$ImageNumber==Image_number]
  
  Dimension = metadata(sce)$dimension
  
  
  if (Dimension == "2D") {
    Temp_size_data = sqrt(sce$Cell_size[sce$ImageNumber==Image_number])
  }
  
  if (Dimension == "3D") {
    Temp_size_data = (sce$Cell_size[sce$ImageNumber==Image_number])^1/3
  }
  
  par(las=1,bty="l")
  plot(Temp_location_data,cex=Temp_size_data/10,pch=21,bg=.cluster_to_color(Temp_cluster_data))
  
}
