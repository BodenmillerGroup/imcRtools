# Helper functions for graphical output
#' @importFrom grDevices colorRamp

.color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("grey","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.99,na.rm = TRUE)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}


.cluster_to_color = function(cluster_vector,Defined_list_cluster = NULL) {
  
  cluster_vector = as.character(cluster_vector)
  List_unique_cluster = unique(cluster_vector)
  List_unique_cluster = List_unique_cluster[order(List_unique_cluster)]
  
  if (!is.null(Defined_list_cluster)) {
    List_unique_cluster = List_unique_cluster
  }
  
  N_clusters = length(unique(cluster_vector))
  
  optimal_palette = suppressWarnings(colorRamp(brewer.pal(N_clusters, "Spectral")))
  optimal_palette = optimal_palette((1:N_clusters)/N_clusters)
  optimal_palette = optimal_palette / 255
  optimal_palette = rgb(optimal_palette)
  
  color_cluster = cluster_vector
  
  for (k in 1:N_clusters) {
    selected_cluster = List_unique_cluster[k]
    color_cluster[cluster_vector==k] = optimal_palette[k]
  }
  return(color_cluster)
}
