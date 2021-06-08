### Helper functions for reading in ImcSegmentationPipeline data
#' @importFrom vroom vroom
.cpout_create_object <- function(path, object_file, image_file, 
                                 object_feature_file, 
                                 intensities, extract_imgid_from,
                                 extract_cellid_from, extract_coords_from,
                                 extract_cellmetadata_from, extract_sampleid_from, 
                                 scale_intensities, extract_scalingfactor_from) {
    
    cur_counts <- vroom(file.path(path, object_file),
                        col_select = c(contains(intensities),
                                       all_of(c(extract_imgid_from,
                                                extract_cellid_from,
                                                extract_coords_from,
                                                extract_cellmetadata_from))))

    scaling_factor <- as.data.frame(vroom(file.path(path, image_file),
                            col_select = all_of(c(extract_imgid_from,
                                                  extract_scalingfactor_from))))
    rownames(scaling_factor) <- as.character(scaling_factor[[extract_imgid_from]])
    
    # Scale counts
    cur_counts <- split(cur_counts, img_id[[extract_imgid_from]])
    scaling_factor <- scaling_factor[names(cur_counts),]
    
    cur_counts <- mapply(function(cells, s_factor) {
        cells %>% mutate(across(contains(intensities), function(x){x * s_factor}))
    }, cur_counts, as.list(scaling_factor[[extract_scalingfactor_from]]), 
    SIMPLIFY = FALSE)
    
    cur_counts <- do.call(rbind, cur_counts)
    
    # Order channels 
    cur_channels <- colnames(cur_counts)[grepl(intensities, colnames(cur_counts))]
    cur_channels_id <- as.numeric(str_extract(cur_channels, "[0-9]{1,3}$"))
    
    cur_channels <- cur_channels[order(cur_channels_id, decreasing = FALSE)]
    
    if (return_as == "spe") {
        object <- SpatialExperiment(assays = list(counts = t(as.matrix(cur_counts[,cur_channels]))),
                                    sample_id = as.character(cur_counts[[extract_imgid_from]]))
        object$ObjectNumber <- cur_counts[[extract_cellid_from]]
    } else {
        object <- SingleCellExperiment(assays = list(counts = t(as.matrix(cur_counts[,cur_channels]))))     
        object$sample_id <- as.character(cur_counts[[extract_imgid_from]])
        object$ObjectNumber <- cur_counts[[extract_cellid_from]]
    }
    
    # Build colData
    if (return_as == "spe") {
        spatialCoords(object) <- matrix(c(cur_counts[[extract_coords_from[1]]],
                                       cur_counts[[extract_coords_from[2]]]),
                                     ncol = 2, byrow = FALSE,
                                     dimnames = list(as.character(object$ObjectNumber),
                                                     c("Pos_X", "Pos_Y")))
        colData(object) <- cbind(colData(object),
                                 cur_counts[,extract_cellmetadata_from])
    } else {
        object$Pos_X <- cur_counts[[extract_coords_from[1]]]
        object$Pos_Y <- cur_counts[[extract_coords_from[2]]]
        colData(object) <- cbind(colData(object),
                                 cur_counts[,extract_cellmetadata_from])
    }
    
    # Set correct metal names
    cur_channels <- vroom(file.path(path, object_feature_file),
                        col_select = c(channel, channel_id))
    cur_channels <- unique(cur_channels)
    cur_channels <- cur_channels[grepl("[a-zA-Z]{1,2}[0-9]{2,3}", cur_channels$channel_id),]
    cur_channels <- cur_channels[order(cur_channels$channel, decreasing = FALSE),]
    
    rownames(object) <- cur_channels$channel_id
}
    
    
.cpout_add_image_metadata <-
    