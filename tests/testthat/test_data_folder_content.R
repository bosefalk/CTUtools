context("Testing data folder content checks")
library(testthat)
library(openssl)

old_wd <- getwd()


# Setup test environment --------------------------------------------------

# Set up test environment in a fresh temp directory, where the files in data/ matches doc/data_folder_content.csv,
# as a baseline we can adjust
baseloc <- tempfile()
dir.create(baseloc)
setwd(baseloc)
dir.create("data")
dir.create("doc")

# Create three data files in the data subdirectory
write.csv(
  data.frame(colA = c("A", "B"), 
             colB = c(13, 90)),
  file = "data/fileA.csv",
  row.names = FALSE
          )
write.table(
  data.frame(colC = c("C", "D"), 
             colD = c(20, 10)),
  file = "data/fileB.txt",
  row.names = FALSE
)
write.csv(
  data.frame(colDate = c(as.Date("2018-01-02"), as.Date("2019-02-02")), 
             colFactor = as.factor(c("factorA", "factorB"))),
  file = "data/fileC.csv",
  row.names = FALSE
)

# Create a data_folder_content matching these files
datafiles <- data.frame(path = list.files("data", recursive = TRUE), stringsAsFactors = FALSE)
for (i in 1:nrow(datafiles)) {
  datafiles$md5[i] <- as.character(openssl::md5(file(paste0("data/", datafiles$path[i]))))
  
}
write.csv(
  datafiles,
  file = "doc/data_folder_content.csv", 
  row.names = FALSE
)

# Function to make a copy of the clean base directory and move the working directory to it
copy_baseloc <- function() {
  tmp_testdir <- tempfile()
  dir.create(tmp_testdir)
  file.copy(paste0(baseloc, "//"), tmp_testdir, recursive = TRUE)
  setwd(tmp_testdir)
}


# Run tests ---------------------------------------------------------------

test_that("Correctly identify missing, changed and renamed files in updated doc/data_folder_content.csv", {
  
  # csv and data files match, no output expected
  copy_baseloc()
  expect_silent(datafolder_check())
  
  
  # Remove one file from data folder
  copy_baseloc()
  file.remove("data/fileB.txt")
  expect_error(datafolder_check(), "data folder file list doesn't match doc/data_folder_content.csv")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "Missing or outdated files in data:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "New files:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "fileB.txt")
  
  
  # Change content of a file
  # Change one cell in fileA.csv, save it in a temp location and update data_folder_content.csv with the md5 from 
  # this new temporary fileA.csv
  copy_baseloc()
  fileA <- read.csv("data/fileA.csv", stringsAsFactors = FALSE)
  fileA$colB[1] <- 15
  tmp_fileA <- tempfile(fileext = ".csv")
  write.csv(fileA, tmp_fileA, row.names = FALSE)
  fileA_newmd5 <- as.character(openssl::md5(file(tmp_fileA)))
  data_folder_content <- read.csv("doc/data_folder_content.csv", stringsAsFactors = FALSE)
  data_folder_content$md5[data_folder_content$path == "fileA.csv"] <- fileA_newmd5
  write.csv(
    data_folder_content,
    file = "doc/data_folder_content.csv", 
    row.names = FALSE
  )
  expect_error(datafolder_check(), "data folder file list doesn't match doc/data_folder_content.csv")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "Missing or outdated files in data:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "Changed files:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "fileA.csv")
  
  
  # Renamed files, rename fileC.csv to fileX.csv in data_folder_content
  copy_baseloc()
  data_folder_content <- read.csv("doc/data_folder_content.csv", stringsAsFactors = FALSE)
  data_folder_content$path[3] <- "fileX.csv"
  write.csv(
    data_folder_content,
    file = "doc/data_folder_content.csv", 
    row.names = FALSE
  )
  expect_error(datafolder_check(), "data folder file list doesn't match doc/data_folder_content.csv")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "Missing or outdated files in data:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "Renamed files:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "fileX.csv")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "fileC.csv")
  


  # One missing file, one changed file and one renamed file:
  copy_baseloc()
  
  file.remove("data/fileB.txt")
  
  fileA <- read.csv("data/fileA.csv", stringsAsFactors = FALSE)
  fileA$colB[1] <- 15
  tmp_fileA <- tempfile(fileext = ".csv")
  write.csv(fileA, tmp_fileA, row.names = FALSE)
  fileA_newmd5 <- as.character(openssl::md5(file(tmp_fileA)))
  data_folder_content <- read.csv("doc/data_folder_content.csv", stringsAsFactors = FALSE)
  data_folder_content$md5[data_folder_content$path == "fileA.csv"] <- fileA_newmd5
  
  data_folder_content$path[3] <- "fileX.csv"
  
  write.csv(
    data_folder_content,
    file = "doc/data_folder_content.csv", 
    row.names = FALSE
  )
  
  expect_error(datafolder_check(), "data folder file list doesn't match doc/data_folder_content.csv")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "Missing or outdated files in data:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "New files:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "fileB.txt")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "Changed files:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "fileA.csv")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "Renamed files:")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "fileX.csv")
  expect_message(suppressWarnings(datafolder_check(stop_on_error = FALSE)), "fileC.csv")
  

  
})


# Return from test environment --------------------------------------------

# Return to the original working directory after the tests have been carried out
setwd(old_wd)
