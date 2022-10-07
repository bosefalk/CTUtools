library(plyr)
library(dplyr)
library(openssl)
library(readxl)

xlsx_file <- "~/Projects/gcp_interlaboratory/data/NMDP DKMS PBSC Graft Composition Study Flow Results_2022 06 24.xlsx"

#' Load FCSExpress output excel
#'
#' Reads the excel sheet output by FCSExpress, optionally calculates the %CD45+ of each population 
#' in addition to the # Events and %Parent already present, and outputs in long format
#'
#' @param xlsx_file file
#' @param CD45 TRUE/FALSE
#' @param format "wide" or "long"
#'
#' @return A data.frame with columns Column and Value as string
#' @seealso \code{\link{HLA_C_class_data}}
#'
#' @examples
#gcp_read_fcsexpress_xlsx <- function(xlsx_file, CD45 = TRUE) {
  
  raw_file_colnames <- suppressMessages(read_excel(xlsx_file, n_max = 0))
  raw_file <- suppressMessages(read_xlsx(xlsx_file, col_types = c(rep("text", 7), rep("numeric", length(raw_file_colnames) - 7))))
  # TODO checks
                                         
                                         
  # There is a naming error in the excel file - two columns with the same name. Following the convention of all other populations
  # and based on the number in the columns, the first one is the Events and the second one %Parent. read_excel automatically appends
  # a numbner to unique columns so they end up ...161 and ...162
  if (all(c("(Non-Classical) Lin- CD11b+ CD33+ CD14- CD16+ Events...161", "(Non-Classical) Lin- CD11b+ CD33+ CD14- CD16+ Events...162") %in%
      colnames(raw_file))) {
    raw_file <- raw_file %>% rename(`(Non-Classical) Lin- CD11b+ CD33+ CD14- CD16+ Events` = `(Non-Classical) Lin- CD11b+ CD33+ CD14- CD16+ Events...161`,
                                    `(Non-Classical) Lin- CD11b+ CD33+ CD14- CD16+ %Parent` = `(Non-Classical) Lin- CD11b+ CD33+ CD14- CD16+ Events...161`)
  }
  
  if (CD45 == TRUE) {
    # 'Events' are now %CD45
    pop_CD45 <- raw_file %>% mutate_at(vars(contains("Events"), -`CD45+ Events`), 
                                       funs(as.numeric(.) / as.numeric(`CD45+ Events`) * 100))
    pop_CD45$`CD45+ Events` <- 100
    
    # Rename Events to %CD45 and add to existing data.frame
    pop_CD45 <- pop_CD45 %>% select(-contains("%Parent")) %>% 
      rename_with(.fn = function(x) gsub("Events", "%CD45", x, fixed = TRUE))
    
    raw_file <- raw_file %>% left_join(pop_CD45)
    
    #"Lin- CD34+ CD38- CD45RA- CD90+ CD117+ ABC", "Lin- CD34+ CD38- CD45RA- CD90+ CD49f+ CD117+ ABC", "Lin- CD34+ CD38- CD45RA- CD90-
    #CD49f+ CD117+ ABC"
  }
  

  
#}