# Helper function for the patch detection method
#' @importFrom sf st_multipoint st_cast st_sfc st_distance
#' @importFrom dplyr as_tibble filter sym nest_by summarize
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
            cells_sfc <- st_cast(st_sfc(cells), "POINT")
            
            if (sum(!is.na(colData(cur_obj)[[name]])) == 0) {
                return(cur_obj)
            }
            
            data <- polygon <- NULL
            
            cur_out <- colData(cur_obj) %>% as_tibble %>%
                filter(!is.na(!!sym(name))) %>% 
                nest_by(!!sym(name)) %>%
                summarize(
                    polygon = list(.polygon_function(x = data, 
                                                    coords = coords, 
                                                    convex = convex)),
                    cells = list(.milieu_function(x = polygon,
                                            distance = expand_by,
                                            cells = cells_sfc)))
            
            # Find cells that are not unique in extended patches
            cur_cells <- do.call(c, cur_out$cells)
            cur_cells <- cur_cells[duplicated(cur_cells)]
            cur_cells <- cur_cells[!is.na(cur_cells)]
            
            if (length(cur_cells) > 0) {
                cur_dists <- mapply(function(y, patch_name){
                    if (is.na(y)) {return(NULL)}
                    cur_mat <- st_distance(cells_sfc[cur_cells], y)
                    colnames(cur_mat) <- patch_name
                    return(cur_mat)
                }, cur_out$polygon, cur_out$patch_id, SIMPLIFY = FALSE)
                cur_dists <- do.call("cbind", cur_dists)
                cur_patch_id <- apply(cur_dists, 1, function(y){
                    return(colnames(cur_dists)[which.min(y)])
                })
            }
            
            cur_patch <- colData(cur_obj)[[name]]
            for (i in seq_len(nrow(cur_out))) {
                if (all(!is.na(cur_out$cells[[i]]))) {
                    cur_patch[cur_out$cells[[i]]] <- cur_out$patch_id[i]
                }
            }
            
            if (length(cur_cells) > 0) {
                cur_patch[cur_cells] <- cur_patch_id
            }
            
            cur_obj[[name]] <- cur_patch
            
            return(cur_obj)
            
        }, BPPARAM = BPPARAM)
    
    return(do.call("cbind", cur_out))
    
}
#' @importFrom concaveman concaveman
#' @importFrom grDevices chull 
#' @importFrom sf st_polygon
.polygon_function <- function(x, coords, convex){
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
        
        return(polygon)
    } else {
        cur_coords <- as.matrix(cbind(x[[coords[1]]], x[[coords[2]]]))
        hull <- data.frame(concaveman(cur_coords, concavity = 1))
        
        polygon <- st_polygon(list(as.matrix(hull)))
        
        return(polygon)
    }
}

#' @importFrom sf st_buffer st_sfc st_intersects
.milieu_function <- function(x, distance, cells){
    
    if (is.na(x)) {
        return(NA)
    }
    
    polygon_buff <- st_buffer(x, distance)
    polygon_buff_sfc <- st_sfc(polygon_buff)
        
    intersect_cells <- st_intersects(polygon_buff_sfc, cells)
        
    return(intersect_cells[[1]])
}
