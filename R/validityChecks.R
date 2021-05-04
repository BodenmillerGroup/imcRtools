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
.validSCEtoTXTinput <- function(txt_list, metadata_cols, verbose){
    
    # Check if input is a named list
    if (!is.list(txt_list)) {
        stop("'txt_list' needs to be a list.")
    }
    
    if (is.null(names(txt_list))) {
        stop("Each entry in 'txt_list' needs to be named.")
    }
    
    # Check if names are all of the format Mt123 where Mt is the metal name
    # and 123 is the mass
    cur_names <- names(txt_list)
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
    cur_channels <- str_extract(colnames(txt_list[[1]]), "[A-Za-z]{1,2}[0-9]{2,3}")
    cur_channels <- cur_channels[!is.na(cur_channels)]
    
    # Verbose option will print possible missmatched between acquired and open
    # channels
    spot_not_ac <- cur_names[!(cur_names %in% cur_channels)]
    ac_not_spot <- cur_channels[!(cur_channels %in% cur_names)]
    if (verbose) {
        message("Spotted channels: ", paste(cur_names, collapse = ", "))
        message("")
        message("Acquired channels: ", paste(cur_channels, collapse = ", "))
        message("")
        message("Channels spotted but not acquired: ", 
                paste(spot_not_ac, collapse = ", "))
        message("")
        message("Channels acquired but not spotted: ",
                paste(ac_not_spot, collapse = ", "))
    }
    
    if (!all(cur_names %in% cur_channels)) {
        stop("Not all spotted channels were acquired.")
    }
    
}

# Check inputs for SCE from TXT function
#' @importFrom SummarizedExperiment colData
.validSpotHeatmapInput <- function(object, spot_id, channel_id, assay_type, log, 
                                   threshold, order_metals){
    
    # Check sce object
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
