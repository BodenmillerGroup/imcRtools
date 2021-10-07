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
            
            
        }, BPPARAM = BPPARAM)
    
}
