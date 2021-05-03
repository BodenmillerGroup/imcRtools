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
