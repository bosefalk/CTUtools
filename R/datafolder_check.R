# For git, we don't want to sync the actual data to gitlab. Instead, the data folder should be manually copied wherever it's safe,
# and there should be a list in doc/data_folder_content.csv with a list of the files expected in data/, and their md5 checksums. 
# This .csv is synced to git, so we can check against it at the start of every analysis script that the actual files in the data 
# folder on this machine matches what the rest of the code expected to be there.
library(openssl)
library(dplyr)

datafolder_check <- function(stop_on_error = TRUE) {
  
  # Check workspace is as expected
  if (any(!c("data", "doc") %in% list.files())) {stop("Missing data and / or doc folder")}
  
  if (!file.exists("doc/data_folder_content.csv")) {
    warning("doc/data_folder_content.csv does not exist so no checking is done, run datafiles_update() to generate list from current data folder")
    return()
  }
  
  if (length(list.files("data")) == 0) {
    warning("data folder is empty, exiting datafiles_check()")
    return()
  }
  
  # Check data_folder_content.csv is as expected
  data_folder_content <- read.csv("doc/data_folder_content.csv", stringsAsFactors = FALSE)
  
  if(!"path" %in% colnames(data_folder_content) | !"md5" %in% colnames(data_folder_content)) {
    stop("path or md5 column missing from doc/data_folder_content.csv")
  }
  
  if (class(data_folder_content$path) != "character" | class(data_folder_content$md5) != "character")
    stop("path or md5 column in doc/data_folder_content.csv not read as character strings")
  
  
  
  # datafiles are the actual files in the data/ folder, data_folder_content is from the doc .csv file
  datafiles <- data.frame(path = list.files("data", recursive = TRUE), stringsAsFactors = FALSE)
  for (i in 1:nrow(datafiles)) {
    datafiles$md5[i] <- as.character(openssl::md5(file(paste0("data/", datafiles$path[i]))))
    
  }
  
  # First look at what appears in the csv but not the data folder
  changes <- dplyr::anti_join(data_folder_content, datafiles, by = c("path", "md5"))

  # Files where the filename doesn't exist are considered new, where the filename exists but md5 is different
  # are changed, and where md5 is same but filename is new are considered renamed
  new_files <- changes %>% filter(!path %in% datafiles$path & !md5 %in% datafiles$md5)
  changed_files <- changes %>% filter(path %in% datafiles$path & !md5 %in% datafiles$md5)
  changed_files <- left_join(changed_files, datafiles, by = "path", suffix = c(".new", ".old"))
  renamed_files <- changes %>% filter(!path %in% datafiles$path & md5 %in% datafiles$md5)
  renamed_files <- left_join(renamed_files, datafiles, by = "md5", suffix = c(".new", ".old"))
  
  # These are files which appear in the data folder but are not listed in data_folder_content
  data_new_files <- dplyr::anti_join(datafiles, data_folder_content, by = c("path", "md5"))
  # Renamed and changed files will already be listed by changed_files and renamed_files above, so only interested in
  # new files which appear in the data folder
  data_new_files <- data_new_files %>% filter(!path %in% data_folder_content$path & !md5 %in% data_folder_content$md5)
  
  # Construct error message to be passed to message()
  error_msg <- function(nf = new_files, cf = changed_files, rf = renamed_files) {
    
    error_string <- "Missing or outdated files in data:\nIf data folder is accurate run CTUtools::datafolder_update()\n"
    
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
  
  if (nrow(data_new_files) + nrow(changes) > 15) {
    error_file <- tempfile()
  }
  
  if (nrow(data_new_files) > 0) {
    dnf_string <- "New files appeared in data folder, run CTUtools::datafolder_update()\n"
    dnf_string <- paste0(dnf_string, paste0(capture.output(data_new_files), collapse = "\n"), "\n")
    
    if (exists("error_file")) {
      sink(error_file)
      cat(dnf_string)
      sink()
      
      if (nrow(changes) == 0) { # if there are only added files the error block which displays the file won't trigger
        file.show(error_file, title = "datafiles_check() messages")
      }
      
    } else {
    message(dnf_string)
    }
  }
  
  # If there are discepancies between the csv and actual data file, send a message with the details and either throw an error or 
  # a warning depending on the stop_on_error flag
  # If the number of discrepancies is too large (>=15) instead of cluttering the console the details are written to a temp file which 
  # is displayed
  if (nrow(changes) > 0) {

    error_string <- error_msg()
    
    if (exists("error_file")) { # Too long to display as a message, open up a file with details of what is missing
      
      sink(error_file, append = TRUE)
      cat(error_string)
      sink()
      file.show(error_file, title = "datafiles_check() messages")
      
      } else {
      message(error_string)
    }

    if (stop_on_error == TRUE) {
      stop("data folder file list doesn't match doc/data_folder_content.csv")
    } else {
      warning("data folder file list doesn't match doc/data_folder_content.csv")
    }
  }
  
  return() # Everything matches
}


datafolder_update <- function() {
  
  # Check workspace is as expected
  if (any(!c("data", "doc") %in% list.files())) {stop("Missing data and / or doc folder")}
  
  if (length(list.files("data")) == 0) {
    warning("data folder is empty, exiting datafiles_update()")
    return()
  }
  
  
  
  
  datafiles <- data.frame(path = list.files("data", recursive = TRUE), stringsAsFactors = FALSE)
  for (i in 1:nrow(datafiles)) {
    datafiles$md5[i] <- as.character(openssl::md5(file(paste0("data/", datafiles$path[i]))))
    
  }
  write.csv(datafiles, file = "doc/data_folder_content.csv", row.names = FALSE)
}