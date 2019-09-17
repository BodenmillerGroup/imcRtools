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