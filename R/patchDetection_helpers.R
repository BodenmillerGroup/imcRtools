# Helper function for the patch detection method
.expand_patch <- function(object, 
                          name,
                          expand_by,
                          coords,
                          convex,
                          img_id,
                          BPPARAM){
    
    cur_out <- bplapply(
        unique(colData(object)[[img_id]]),
        function(x){
            
            cur_obj <- object[,as.character(colData(object)[[img_id]]) == x]
            
            if (is(cur_obj, "SpatialExperiment")) {
                cur_coords <- spatialCoords(cur_obj)[,coords]
            } else {
                cur_coords <- colData(cur_obj)[,coords]
            }
            
            cells <- st_multipoint(as.matrix(cur_coords))
            cells_sfc = st_cast(st_sfc(cells), "POINT")
            
            cur_out <- colData(cur_obj) %>% as_tibble %>%
                filter(!is.na(!!sym(name))) %>% 
                nest_by(!!sym(name)) %>%
                summarize(polygon = list(hull_function(x = data, 
                                                  coords = coords, 
                                                  convex = convex)))
            
            
        }, BPPARAM = BPPARAM)
    
}

hull_function <- function(x, coords, convex){
    
    if (convex) {
        hull <- chull(x = x[[coords[1]]], y = x[[coords[2]]])
        
        # cells that build the border of a patch
        border_cells = x[hull,]
        coordinates = as.matrix(border_cells[,coords])
        coordinates <- rbind(coordinates, coordinates[1,])

        return(st_polygon(list(coordinates)))
    } else {
        cur_coords <- as.matrix(cbind(x[[coords[1]]], x[[coords[2]]]))
        hull <- data.frame(concaveman(cur_coords, concavity = 1))

        return(st_polygon(list(as.matrix(hull))))
    }
    
}
