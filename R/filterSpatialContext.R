filterSpatialContext <- function(object,
                               entry = "spatial_context",
                               img_id = "sample_id",
                               n_cells_threshold = NULL,
                               n_samples_threshold = NULL,
                               name = "spatial_context_filtered"
                               ){
  
  #.valid.filterSpatialContext.input(object, entry, img_id, n_cells, n_samples) TO DO
  
  data <- colData(object)[,colnames(colData(object)) %in% c(entry,img_id)] %>% 
  table() %>% as.data.frame
  Freq <- as.name("Freq")
  n <- as.name("n")
  
  anno <- data.frame(spatial_context = unique(data[,entry]), 
                     n_cells = data %>% group_by_at(entry) %>% 
                       summarise(sum = sum(Freq)) %>% pull(sum), 
                     n_samples = data %>% group_by_at(entry) %>% 
                       filter(Freq != 0) %>% count() %>% pull(n) %>% 
                       as.integer()
                     )
  
  #option for filtering just one level - TODO (if etc.)
  selected <- anno %>% 
    filter(n_cells >= n_cells_threshold & n_samples >= n_samples_threshold) %>% 
    pull(spatial_context) %>% unfactor()       
  
  out_dat <- ifelse(colData(object)[[entry]] %in% selected, 
                    colData(object)[[entry]], NA)
  
  colData(object)[[name]] <- out_dat
  return(object)
}
