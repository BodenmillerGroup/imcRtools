#' @rdname plot_embedding
#' @title Plot a low-dimensional embedding of single-cell data
#'
#' @description Plot a low-dimensional embedding of single-cell data that has been previously computed and is stored in a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param selected_embedding type of graph to compute. Has to be chosen among the following values : "KNN","FRNN" or "Delaunay"
#' @param overlay_cluster boolean value. If set to TRUE, will color the points according to their clustering 
#' @return Return a colored plot of a low-dimensional embedding 
#' @examples 
#' sce = plot_embedding(sce,selected_embedding = "UMAP",overlay_cluster=T)
#'
#'@import SingleCellExperiment


plot_embedding = function(sce,selected_embedding = "UMAP",overlay_cluster=F) {
  if (!selected_embedding%in%reducedDimNames(sce)) {
    stop("The selected embedding has not been computed. Please compute it first or check that you have chosen the correct embedding \n")
  }
  
  Embedding = reducedDim(sce)
  if (class(Embedding)=="list") {
    Embedding = Embedding[["selected_embedding"]]
  }
  par(las=1,bty="l")
  
  if (is.null(colLabels(sce))) {
    plot(Embedding,pch=21,bg="red3", xlab="Dimension 1",ylab="Dimension 2",cex=0.7)
  }
  
  if (!is.null(colLabels(sce))) {
    plot(Embedding,pch=21,bg=.cluster_to_color(colLabels(sce)), xlab="Dimension 1",ylab="Dimension 2",cex=0.7)
  }
}
