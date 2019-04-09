# For git, we don't want to sync the actual data to gitlab. Instead, the data folder should be manually copied wherever it's safe,
# and there should be a list in doc/data_folder_content.csv with a list of the files expected in data/, and their md5 checksums. 
# This .csv is synced to git, so we can check against it at the start of every analysis script that the actual files in the data 
# folder on this machine matches what the rest of the code expected to be there.
library(openssl)
library(dplyr)

datafolder_check <- function(stop_on_error = TRUE) {
  # datafiles are the actual files in the data/ folder, data_folder_content is from the doc .csv file
  datafiles <- data.frame(path = list.files("data", recursive = TRUE), stringsAsFactors = FALSE)
  for (i in 1:nrow(datafiles)) {
    datafiles$md5[i] <- as.character(openssl::md5(file(paste0("data/", datafiles$path[i]))))
    
  }
  
  data_folder_content <- read.csv("doc/data_folder_content.csv", stringsAsFactors = FALSE)
  
  # First look at what appears in the csv but not the data folder
  changes <- dplyr::anti_join(data_folder_content, datafiles, by = c("path", "md5"))

  # Files where the filename doesn't exist are considered new, where the filename exists but md5 is different
  # are changed, and where md5 is same but filename is new are considered renamed
  new_files <- changes %>% filter(!path %in% datafiles$path & !md5 %in% datafiles$md5)
  changed_files <- changes %>% filter(path %in% datafiles$path & !md5 %in% datafiles$md5)
  changed_files <- left_join(changed_files, datafiles, by = "path", suffix = c(".new", ".old"))
  renamed_files <- changes %>% filter(!path %in% datafiles$path & md5 %in% datafiles$md5)
  renamed_files <- left_join(renamed_files, datafiles, by = "md5", suffix = c(".new", ".old"))
  
  # Construct error message to be passed to message()
  error_msg <- function(nf = new_files, cf = changed_files, rf = renamed_files) {
    
    error_string <- "Missing or outdated files in data:\n"
    
    if (nrow(nf) > 0) {
      error_string <- paste0(error_string, "New files:\n")
      error_string <- paste0(error_string, paste0(capture.output(new_files), collapse = "\n"), "\n")
    }
    
    if (nrow(cf) > 0) {
      error_string <- paste0(error_string, "Changed files:\n")
      error_string <- paste0(error_string, paste0(capture.output(changed_files), collapse = "\n"), "\n")
      
    }
    
    if (nrow(rf) > 0) {
      error_string <- paste0(error_string, "Renamed files:\n")
      error_string <- paste0(error_string, paste0(capture.output(renamed_files), collapse = "\n"), "\n")
    }
    
    return(error_string)
  }
  
  # If there are discepancies between the csv and actual data file, send a message with the details and either throw an error or 
  # a warning depending on the stop_on_error flag
  # If the number of discrepancies is too large (>=15) instead of cluttering the console the details are written to a temp file which 
  # is displayed
  if (nrow(changes) > 0) {

    error_string <- error_msg()
    
    if (nrow(changes) > 15) { # Too long to display as a message, open up a file with details of what is missing
      error_file <- tempfile()
      sink(error_file)
      cat(error_string)
      sink()
      file.show(error_file)
    } else {
      message(error_string)
    }
    
    if (stop_on_error == TRUE) {
      stop("data folder file list doesn't match doc/data_folder_content.csv")
    } else {
      warning("data folder file list doesn't match doc/data_folder_content.csv")
    }
  }
  
  
}