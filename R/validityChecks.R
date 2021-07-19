# Checks if input files exist and if they contain the right entries
.fileChecks <- function(cells, image_meta, relationships, panel){
  
  # Object to save errors
  errors <- c()
  
  # Cells.csv file
  if(!file.exists(cells)){
    errors <- c(errors, paste0("The file containing the cell features does not exist. ", 
                "Please correct the path to the file.\n"))
  }
  
  # Image metadata file
  if(!file.exists(image_meta)){
    errors <- c(errors, paste0("The file containing the image metadata does not exist. ",
                "Please correct the path to the file.\n"))
  }
  
  # Relationship file
  if(!file.exists(relationships)){
    errors <- c(errors, paste0("The file containing the cell-cell relationships does not exist. ",
                "Please correct the path to the file.\n"))
  }
  
  # Panel file
  if(!file.exists(panel)){
    errors <- c(errors, paste0("The file containing the panel does not exist. ",
                "Please correct the path to the file.\n"))
  }
  
  return(errors)
}

# Check function for SCE from TXT function
#' @importFrom stringr str_extract
.valid.readSCEfromTXT.input <- function(txt_list, cur_names, 
                                        metadata_cols, verbose){
    
    # Check if names are all of the format Mt123 where Mt is the metal name
    # and 123 is the mass
    cur_mass <- str_extract(cur_names, "[0-9]{2,3}$")
    cur_names <- cur_names[order(as.numeric(cur_mass))]
    
    if (!all(grepl("^[A-Za-z]{2}[0-9]{2,3}$", cur_names))) {
        stop("Not all names match the pattern (mt)(mass).")
    }
    
    # Check if metadata_cols are in each file
    cur_check <- lapply(txt_list, function(x){
        all(metadata_cols %in% colnames(x))
    })
    
    if (!all(unlist(cur_check))) {
        stop("Not all 'metadata_cols' are present in entries to 'txt_list'")
    }
    
    # Check if spotted channel is also open
    cur_channels <- str_extract(colnames(txt_list[[1]]), 
                                "[A-Za-z]{1,2}[0-9]{2,3}Di")
    cur_channels <- cur_channels[!is.na(cur_channels)]
    cur_channels <- sub("Di", "", cur_channels)
    
    # Verbose option will print possible missmatched between acquired and open
    # channels
    spot_not_ac <- cur_names[!(cur_names %in% cur_channels)]
    ac_not_spot <- cur_channels[!(cur_channels %in% cur_names)]
    if (verbose) {
        cat("Spotted channels: ", paste(cur_names, collapse = ", "))
        cat("\n")
        cat("Acquired channels: ", paste(cur_channels, collapse = ", "))
        cat("\n")
        cat("Channels spotted but not acquired: ", 
                paste(spot_not_ac, collapse = ", "))
        cat("\n")
        cat("Channels acquired but not spotted: ",
                paste(ac_not_spot, collapse = ", "))
    }
    
    if (!all(cur_names %in% cur_channels)) {
        stop("Not all spotted channels were acquired.")
    }
    
}

#' @importFrom SummarizedExperiment colData assayNames
.valid.plotSpotHeatmap.input <- function(object, spot_id, channel_id, 
                                         assay_type, log, 
                                   threshold, order_metals){
    
    # Check sce object
    if (!is(object, "SingleCellExperiment")) {
        stop("'object' needs to be a SingleCellExperiment object.")
    }
    if (!spot_id %in% names(colData(object))) {
        stop("'spot_id' not in 'colData(object)'.")
    }
    if (!channel_id %in% names(rowData(object))) {
        stop("'channel_id' not in 'rowData(object)'.")
    }
    if (!assay_type %in% assayNames(object)) {
        stop("'assay_type' not in 'assayNames(object)'.")
    }
    
    if (!all(is.logical(log)) || length(log) != 1) {
        stop("'log' needs to be logical.")
    }
    
    if (!is.null(threshold) & (length(threshold) != 1 || !all(is.numeric(threshold)))) {
        stop("'threshold' needs to be a single numeric.")
    }
    
    if (!all(is.logical(order_metals)) || length(order_metals) != 1) {
        stop("'order_metals' needs to be logical.")
    }
}

#' @importFrom SummarizedExperiment colData assayNames
.valid.binAcrossPixels.input <- function(object, bin_size, spot_id, assay_type){
    
    # Check sce object
    if (!is(object, "SingleCellExperiment")) {
        stop("'object' needs to be a SingleCellExperiment object.")
    }
    if (!spot_id %in% names(colData(object))) {
        stop("'spot_id' not in 'colData(object)'.")
    }
    if (!assay_type %in% assayNames(object)) {
        stop("'assay_type' not in 'assayNames(object)'.")
    }
    
    if (!all(is.numeric(bin_size)) || length(bin_size) != 1) {
        stop("'bin_size' needs to be single numeric.")
    }
}

.valid.filterPixels.input <- function(object, bc_id, spot_mass, 
                                    minevents, correct_pixels){
    
    # Check sce object
    if (!is(object, "SingleCellExperiment")) {
        stop("'object' needs to be a SingleCellExperiment object.")
    }
    if (!bc_id %in% names(colData(object))) {
        stop("'bc_id' not in 'colData(object)'.")
    }
    if (!spot_mass %in% names(colData(object))) {
        stop("'spot_mass' not in 'colData(object)'.")
    }

    
    if (!all(is.numeric(minevents)) || length(minevents) != 1) {
        stop("'minevents' needs to be single numeric.")
    }
    if (!all(is.logical(correct_pixels)) || length(correct_pixels) != 1) {
        stop("'correct_pixels' needs to be single logical.")
    }
}

.valid.readImagefromTXT.input <- function(path, pattern){
    # Check if input is character
    if (!all(is.character(path)) | length(path) != 1) {
        stop("Please provide a string input indicating a path.")
    }
    
    if (!is.character(pattern) | length(pattern) != 1) {
        stop("'pattern' must be indicated as single character")
    }
    
    if (!dir.exists(path)) {
        stop("The indicated path does not exist")
    }
            
    out <- list.files(path, pattern = pattern, full.names = TRUE)
    
    if (length(out) == 0) {
        stop("The pattern does not match any\n",
             "of the files in the provided directory.")
    }
}


.valid.read_steinbock.input <- function(path, intensities_folder, graphs_folder,
                                        regionprops_folder, cell_id, coords,
                                        panel, name, pattern) {
    if (length(path) != 1 | !is.character(path)) {
        stop("'path' must be a single string.")
    }
    
    if (!dir.exists(path)) {
        stop("'path' doesn't exist.")
    }
    
    if (is.null(intensities_folder)) {
        stop("'intensities_folder' must be specified.")
    }
    
    if (length(intensities_folder) != 1 | !is.character(intensities_folder)) {
        stop("'intensities_folder' must be a single string.")
    }
    
    if (!dir.exists(file.path(path, intensities_folder))) {
        stop("'intensities_folder' doesn't exist.")
    }
    
    if (!is.null(graphs_folder)) {
        
        if (length(graphs_folder) != 1 | !is.character(graphs_folder)) {
            stop("'graphs_folder' must be a single string.")
        }
        
        if (!dir.exists(file.path(path, graphs_folder))) {
            stop("'graphs_folder' doesn't exist.")
        }
        
    }
    
    if (!is.null(regionprops_folder)) {
        
        if (length(regionprops_folder) != 1 | !is.character(regionprops_folder)) {
            stop("'regionprops_folder' must be a single string.")
        }
        
        if (!dir.exists(file.path(path, regionprops_folder))) {
            stop("'regionprops_folder' doesn't exist.")
        }
        
    }
    
    # Check if any files can be read in
    all_files <- list.files(file.path(path, intensities_folder), 
                            pattern =  pattern, full.names = TRUE)
    
    if (length(all_files) == 0) {
        stop("No files were read in.")
    }
    
    # Check cell_id
    cur_file <- vroom(all_files[1],
                      progress = FALSE, 
                      col_types = cols())
    
    if (length(cell_id) != 1 | !is.character(cell_id)) {
        stop("'extract_cellid_from' must be a single string.")
    }
    
    if (!cell_id %in% colnames(cur_file)) {
        stop("'extract_cellid_from' not in intensities files.")
    }
    
    # Check coords    
    if (!is.null(regionprops_folder)) {
        all_files <- list.files(file.path(path, regionprops_folder), 
                                pattern =  pattern, full.names = TRUE)
        cur_file <- vroom(all_files[1],
                          progress = FALSE, 
                          col_types = cols())
        
        if (!all(is.character(coords))) {
            stop("'extract_coords_from' must be characters.")
        }
        
        if (!all(coords %in% colnames(cur_file))) {
            stop("'coords' not in regionprops files.")
        }
    }
    
    # Check panel
    if (!is.null(panel)) {
        if (length(panel) != 1 | !is.character(panel)) {
            stop("'panel_file' must be a single string.")
        }
        
        if (file.exists(file.path(path, panel))) {
            cur_panel <- vroom(file.path(path, panel), 
                               progress = FALSE, 
                               col_types = cols())
        } else if (file.exists(panel)) {
            cur_panel <- vroom(panel, 
                               progress = FALSE, 
                               col_types = cols())
        } 
        
        if (exists("cur_panel")) {
            
            if (length(name) != 1 | !is.character(name)) {
                stop("'extract_names_from' must be a single string.")
            }
            
            if (!name %in% colnames(cur_panel)) {
                stop("'extract_names_from' not in panel file.")
            }
        }
    }
}

#' @importFrom stringr str_count
.valid.read_cpout.input <- function(path, object_file, image_file,
                                    panel_file, graph_file, object_feature_file,
                                    intensities, 
                                    extract_imgid_from, extract_cellid_from, 
                                    extract_coords_from,
                                    extract_cellmetadata_from, extract_imagemetadata_from,
                                    extract_graphimageid_from, extract_graphcellids_from,
                                    extract_metal_from, scale_intensities,
                                    extract_scalingfactor_from){
    
    if (length(path) != 1 | !is.character(path)) {
        stop("'path' must be a single string.")
    }
    
    if (!dir.exists(path)) {
        stop("'path' doesn't exist.")
    }
    
    if (is.null(object_file)) {
        stop("'object_file' must be specified.")
    }
    
    if (is.null(object_feature_file)) {
        stop("'object_feature_file' must be specified.")
    }
    
    if (length(object_file) != 1 | !is.character(object_file)) {
        stop("'object_file' must be a single string.")
    }
    
    if (length(object_feature_file) != 1 | !is.character(object_feature_file)) {
        stop("'object_feature_file' must be a single string.")
    }
    
    if (!file.exists(file.path(path, object_file))) {
        stop("'object_file' doesn't exist.")
    }
    
    if (!file.exists(file.path(path, object_feature_file))) {
        stop("'object_feature_file' doesn't exist.")
    }
    
    if (!is.null(image_file)) {
        
        if (length(image_file) != 1 | !is.character(image_file)) {
            stop("'image_file' must be a single string.")
        }
        
        if (!file.exists(file.path(path, image_file))) {
            stop("'image_file' doesn't exist.")
        }
        
    }
    
    if (!is.null(graph_file)) {
        
        if (length(graph_file) != 1 | !is.character(graph_file)) {
            stop("'graph_file' must be a single string.")
        }
        
        if (!file.exists(file.path(path, graph_file))) {
            stop("'graph_file' doesn't exist.")
        }
        
    }
    
    # Check panel
    if (!is.null(panel_file)) {
        if (length(panel_file) != 1 | !is.character(panel_file)) {
            stop("'panel_file' must be a single string.")
        }
        
        if (file.exists(file.path(path, panel_file))) {
            cur_panel <- vroom(file.path(path, panel_file), 
                               progress = FALSE, 
                               col_types = cols())
        } else if (file.exists(panel_file)) {
            cur_panel <- vroom(panel_file, 
                               progress = FALSE, 
                               col_types = cols())
        } 
        
        if (exists("cur_panel")) {
            
            if (is.null(extract_metal_from)) {
                stop("'extract_metal_from' must be specified.")
            }
            
            if (length(extract_metal_from) != 1 | !is.character(extract_metal_from)) {
                stop("'extract_metal_from' must be a single string.")
            }
            
            if (!extract_metal_from %in% colnames(cur_panel)) {
                stop("'extract_metal_from' not in panel file.")
            }
        }
    }
    
    # Check object files
    cur_file <- vroom(file.path(path, object_file), n_max = 1, col_types = cols())
    
    if (is.null(intensities)) {
        stop("'intensities' must be specified.")
    }
    
    if (length(intensities) != 1 | !is.character(intensities)) {
        stop("'intensities' must be a single string.")
    }
    
    feature_count <- sum(str_count(colnames(cur_file), intensities))
    
    if (feature_count == 0) {
        stop("No intensity features were read in. Please check the 'intensities' parameter.")
    }
    
    cur_features <- cur_file %>% select(contains(intensities))
    cur_channels <- str_extract(colnames(cur_features), "c[0-9]*$")
    cur_channels <- as.numeric(table(cur_channels))
    
    if (any(cur_channels > 1)) {
        stop("Some of the features set via 'intensities' cannot be uniquely accessed.")
    }
    
    if (is.null(extract_imgid_from)) {
        stop("'extract_imgid_from' must be specified.")
    }
    
    if (length(extract_imgid_from) != 1 | !is.character(extract_imgid_from)) {
        stop("'extract_imgid_from' must be a single string.")
    }
    
    if (!extract_imgid_from %in% colnames(cur_file)) {
        stop("'extract_imgid_from' not in 'object_file'.")
    }
    
    if (is.null(extract_cellid_from)) {
        stop("'extract_cellid_from' must be specified.")
    }
    
    if (length(extract_cellid_from) != 1 | !is.character(extract_cellid_from)) {
        stop("'extract_cellid_from' must be a single string.")
    }
    
    if (!extract_cellid_from %in% colnames(cur_file)) {
        stop("'extract_cellid_from' not in 'object_file'.")
    }
    
    if (!is.null(extract_coords_from)) {
        if (!all(extract_coords_from %in% colnames(cur_file))) {
            stop("'extract_coords_from' not in 'object_file'.")
        }
    }
    
    if (!is.null(extract_cellmetadata_from)) {
        if (!all(extract_cellmetadata_from %in% colnames(cur_file))) {
            stop("'extract_cellmetadata_from' not in 'object_file'.")
        }
    }
    
    
    # Check image files
    if (is.null(scale_intensities)) {
        stop("'scale_intensities' needs to be logical.")
    }
    
    if (length(scale_intensities) != 1) {
        stop("'scale_intensities' needs to be of length 1.")
    } else {
        if (!is.logical(scale_intensities)) {
            stop("'scale_intensities' needs to be logical.")
        }
    }
    
    if (scale_intensities & is.null(image_file)) {
        stop("When scaling the summarized object intensities, ",
        "please supply the 'image_file'.")
    }
    
    if (!is.null(image_file)) {
        cur_file <- vroom(file.path(path, image_file), n_max = 1, col_types = cols())
        
        if (!is.null(extract_imagemetadata_from)) {
            if (!all(extract_imagemetadata_from %in% colnames(cur_file))) {
                stop("'extract_imagemetadata_from' not in 'image_file'.")
            }
        }
        
        if (length(extract_scalingfactor_from) != 1 | !is.character(extract_scalingfactor_from)) {
            stop("'extract_scalingfactor_from' must be a single string.")
        }
        
        if (scale_intensities & !extract_scalingfactor_from %in% colnames(cur_file)) {
            stop("'extract_scalingfactor_from' not in 'image_file'.")
        }
        
    }
    
    # Check graph file
    if (!is.null(graph_file)) {
        cur_file <- vroom(file.path(path, graph_file), n_max = 1, col_types = cols())
        
        if (is.null(extract_graphimageid_from)) {
            stop("'extract_graphimageid_from' must be specified.")
        }
        
        if (length(extract_graphimageid_from) != 1 | !is.character(extract_graphimageid_from)) {
            stop("'extract_graphimageid_from' must be a single string.")
        }
        
        if (!extract_graphimageid_from %in% colnames(cur_file)) {
            stop("'extract_graphimageid_from' not in 'graph_file'.")
        }
        
        if (is.null(extract_graphcellids_from)) {
            stop("'extract_graphcellids_from' must be specified.")
        }
        
        if (!all(extract_graphcellids_from %in% colnames(cur_file))) {
            stop("'extract_graphcellids_from' not in 'graph_file'.")
        }
        
    }
    
}

#' @importFrom SpatialExperiment spatialCoordsNames
.valid.buildSpatialGraph.input <- function(object, type, img_id, k, threshold, coords,
                                    name, directed){
    
    if (!is(object, "SingleCellExperiment")) {
        stop("'object' not of type 'SingleCellExperiment'.")
    }
    
    if (length(img_id) != 1 | !is.character(img_id)) {
        stop("'img_id' must be a single string.")
    }
    
    if (!img_id %in% names(colData(object))) {
        stop("'img_id' not in colData(object).")
    }
    
    if (type == "expansion") {
        
        if (length(threshold) != 1) {
            stop("'threshold' must be a single numeric")
        }
        
        if (is.null(threshold)) {
            stop("When constructing a graph via expansion, please specify 'threshold'.")
        }
    }
    
    if (type == "knn") {
        
        if (length(k) != 1) {
            stop("'k' must be a single numeric")
        }
        
        if (is.null(k)) {
            stop("When constructing a graph via nearest neighbour detection, ",
                 "please specify 'k'.")
        }
    }
    
    if (length(coords) != 2 | !all(is.character(coords))) {
        stop("'coords' must be a character vector of length 2.")
    }
    
    if (is(object, "SpatialExperiment")) {
        if (!all(coords %in% spatialCoordsNames(object))) {
            stop("'coords' not in spatialCoords(object).")
        }
    } else {
        if (!all(coords %in% names(colData(object)))) {
            stop("'coords' not in colData(object).")
        }
    }
    
    if (!is.null(name)) {
        if (length(name) != 1 | !is.character(name)) {
            stop("'name' must be a single string.")
        }
    }
    
    if (length(directed) != 1 | !is.logical(directed)) {
        stop("'directed' must be a single logical.")
    }
    
}

#' @importFrom S4Vectors mcols
.valid.plotSpatial.input <- function(object, img_id, coords, node_color_by, 
                                     node_shape_by, node_size_by, edge_color_by,
                                     edge_width_by, draw_edges, directed, arrow, 
                                     end_cap, colPairName, ncols, nrows, scales){
    
    if (!is(object, "SingleCellExperiment")) {
        stop("'object' not of type 'SingleCellExperiment'.")
    }
    
    if (length(img_id) != 1 | !is.character(img_id)) {
        stop("'img_id' must be a single string.")
    }
    
    if (!img_id %in% names(colData(object))) {
        stop("'img_id' not in colData(object).")
    }
    
    if (length(coords) != 2 | !all(is.character(coords))) {
        stop("'coords' must be a character vector of length 2.")
    }
    
    if (is(object, "SpatialExperiment")) {
        if (!all(coords %in% spatialCoordsNames(object))) {
            stop("'coords' not in spatialCoords(object).")
        }
    } else {
        if (!all(coords %in% names(colData(object)))) {
            stop("'coords' not in colData(object).")
        }
    }
    
    if (!is.null(node_color_by) && (length(node_color_by) != 1 | 
        !is.character(node_color_by))) {
        stop("'node_color_by' must be a single string.")
    }
    
    if (!is.null(node_color_by) &&
        (!node_color_by %in% names(colData(object)))) {
        stop("'node_color_by' not in colData(object).")
    }
    
    if (!is.null(node_shape_by) && 
        (length(node_shape_by) != 1 | !is.character(node_shape_by))) {
        stop("'node_shape_by' must be a single string.")
    }
    
    if (!is.null(node_shape_by) && 
        (!node_shape_by %in% names(colData(object)))) {
        stop("'node_shape_by' not in colData(object).")
    }
    
    if (!is.null(node_size_by) &&
        (length(node_size_by) != 1 | !is.character(node_size_by))) {
        stop("'node_size_by' must be a single string.")
    }
    
    if (!is.null(node_size_by) &&
        (!node_size_by %in% names(colData(object)))) {
        stop("'node_size_by' not in colData(object).")
    }
    
    if (length(draw_edges) != 1 | !is.logical(draw_edges)) {
        stop("'draw_edges' must be a single logical")
    }
    
    if (draw_edges) {
        if (is.null(colPairName)) {
            stop("Please specify the name of the", 
                 " column pairing via 'colPairName'.")
        }
        
        if (length(colPairName) != 1 | !is.character(colPairName)) {
            stop("'colPairName' must be a single string.")
        }
        
        if (!colPairName %in% colPairNames(object)) {
            stop("No column pairing with name ", colPairName, " found.")
        }
        
        if (!is.null(edge_color_by) && 
            (length(edge_color_by) != 1 | !is.character(edge_color_by))) {
            stop("'edge_color_by' must be a single string.")
        }
        
        if (!is.null(edge_color_by) &&
            (!edge_color_by %in% names(colData(object)) &&
            !edge_color_by %in% names(mcols(colPair(object, colPairName))))) {
            stop("'edge_color_by' not in 'colData(object)'", 
                 " or in 'mcols(colPair(object, colPairName))'.")
        }
        
        if (!is.null(edge_width_by) && 
            (length(edge_width_by) != 1 | !is.character(edge_width_by))) {
            stop("'edge_width_by' must be a single string.")
        }
        
        if (!is.null(edge_width_by) && 
            (!edge_width_by %in% names(colData(object)) &&
            !edge_width_by %in% names(mcols(colPair(object, colPairName))))) {
            stop("'edge_width_by' not in 'colData(object)'", 
                 " or in 'mcols(colPair(object, colPairName))'.")
        }
        
        if (length(directed) != 1 | !is.logical(directed)) {
            stop("'directed' must be a single logical")
        }
        
        if (!is.null(arrow) && !is(arrow, "arrow")) {
            stop("'arrow' must be of class grid::arrow.")
        }
        
        if (!is.null(end_cap) && !is(end_cap, "geometry")) {
            stop("'end_cap' must be of type 'geometry'.")
        }
    }
    
    if (!is.null(ncols) && (length(ncols) != 1 | !is.numeric(ncols))) {
        stop("'ncols' must be a single numeric")
    }
    
    if (!is.null(nrows) && (length(nrows) != 1 | !is.numeric(nrows))) {
        stop("'nrows' must be a single numeric")
    }
    
    if (!scales %in% c("fixed", "free_x", "free_y", "free")) {
        stop("'scales' should be one of 'fixed', 'free_x', 'free_y', 'free'.")
    }
    
    
}
