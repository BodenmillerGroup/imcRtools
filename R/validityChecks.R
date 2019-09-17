# Checks if input files exist and if they contain the right entries
.fileChecks <- function(cells, image_meta, relationships, panel){
  
  # Object to save errors
  errors <- c()
  
  # Cells.csv file
  if(!file.exists(cells)){
    errors <- c(errors, "The file containing the cell features does not exist. \n 
                Please correct the path to the file.")
  }
  
  # Image metadata file
  if(!file.exists(image_meta)){
    errors <- c(errors, "The file containing the image metadata does not exist. \n 
                Please correct the path to the file.")
  }
  
  # Relationship file
  if(!file.exists(relationships)){
    errors <- c(errors, "The file containing the cell-cell relationships does not exist. \n 
                Please correct the path to the file.")
  }
  
  # Panel file
  if(!file.exists(panel)){
    errors <- c(errors, "The file containing the panel does not exist. \n 
                Please correct the path to the file.")
  }
  
  return(errors)
}