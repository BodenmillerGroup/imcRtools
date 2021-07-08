# function to generate a data frame with the count of cells of a specified type (e.g. clustered celltypes) in the neighborhood of each cell.
# ... to be written

summarizeEnvironment <- function(object,
                                 graphName,
                                 group,
                                 name = NULL){

  if (! graphName %in% colPairNames(object)) {
    stop("graphName is not a graph in the current object")
  }

  if (! group %in% colnames(colData(object))) {
    stop("group is not a valid colData entry in the current object")
  }

  cur_graph <- graph_from_edgelist(as.matrix(colPair(object,graphName)))

  cur_dat <- as_data_frame(cur_graph)

  cur_out <- cur_dat %>%
    mutate(celltype =   colData(object)[[group]][cur_dat$to]) %>%
    group_by(from,.drop=FALSE) %>%
    count(celltype) %>%
    pivot_wider(names_from = celltype,values_from = n, values_fill = 0) %>%
    ungroup() %>%
    select(unique(colData(object)[[group]]))

  # add the dataframe to the object

  name <- ifelse(is.null(name), "summarizedEnvironment", name)

  colData(object)[[name]] <- cur_out

  return(object)

}
