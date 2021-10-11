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
            
            if (sum(!is.na(colData(cur_obj)[[name]])) == 0) {
                return(cur_obj)
            }
            
            cur_out <- colData(cur_obj) %>% as_tibble %>%
                filter(!is.na(!!sym(name))) %>% 
                nest_by(!!sym(name)) %>%
                summarize(cells = list(milieu_function(x = data, 
                                                  coords = coords, 
                                                  convex = convex,
                                                  distance = expand_by,
                                                  cells = cells_sfc)))
            
            cur_patch <- colData(cur_obj)[[name]]
            for (i in nrow(cur_out)) {
                if (!all(is.na(cur_out$cells[[i]]))) {
                    cur_patch[cur_out$cells[[i]]] <- cur_out$patch_id[i]
                }
            }
            
            cur_obj[[name]] <- cur_patch
            
            return(cur_obj)
            
        }, BPPARAM = BPPARAM)
    
    return(do.call("cbind", cur_out))
    
}

milieu_function <- function(x, coords, convex, distance, cells){
    
    if (nrow(x) <= 2) {
        return(NA)
    }
    
    if (convex) {
        hull <- chull(x = x[[coords[1]]], y = x[[coords[2]]])
        
        # cells that build the border of a patch
        border_cells = x[hull,]
        coordinates = as.matrix(border_cells[,coords])
        coordinates <- rbind(coordinates, coordinates[1,])

        polygon <- st_polygon(list(coordinates))
        
        polygon_buff <- st_buffer(polygon, distance)
        polygon_buff_sfc <- st_sfc(polygon_buff)
        
        intersect_cells <- st_intersects(polygon_buff_sfc, cells)
        
        return(intersect_cells[[1]])
    } else {
        cur_coords <- as.matrix(cbind(x[[coords[1]]], x[[coords[2]]]))
        hull <- data.frame(concaveman(cur_coords, concavity = 1))

        polygon <- st_polygon(list(as.matrix(hull)))
        
        polygon_buff <- st_buffer(polygon, distance)
        polygon_buff_sfc <- st_sfc(polygon_buff)
        
        intersect_cells <- st_intersects(polygon_buff_sfc, cells)
        
        return(intersect_cells[[1]])
    }
    
}
