#' Reduce HLA allele string to specified number of fields
#'
#' Reduce a long HLA allele string (i.e. "01:02:01:03g") down to fewer fields,
#' default being 2 fields. NMDP codes are preserved, trailing "G" are removed, NA or empty string input generates NA output,
#' anything else is output as-is, cast as string. 
#' No input or output validation is done on the data to ensure it meets expected HLA allele string formats. 
#'
#' @param allele String, including ":" but not initial "A*", so for example "01:02:01:03"
#' @param fields Number of fields to keep, default = 2 which would return "01:02"
#'
#' @return character vector with shortened fields
#'
#' @examples
#' dat <- data.frame(HLA_A1 = c("01:02:01:03", "03:ADJRE", "01:02:01:01G"), stringsAsFactors = FALSE)
#' shorten_allele(dat$HLA_A1)
#' shorten_allele(dat$HLA_A1, fields = 3)
#' dat2 <- data.frame(HLA_A1 = c("01:02:01:03", "03:ADJRE", NA, "some_random_string", 2, ""), stringsAsFactors = FALSE)
#' shorten_allele(dat2$HLA_A1)
shorten_allele <- function(allele, fields = 2) {


  if (!fields %in% c(1,2,3,4)) {
    stop("shortern_allele: fields can only be set to 1,2,3 or 4")
  }



  # Remove trailing g if present
  allele <- gsub("G$", "", allele)

  split_allele <- strsplit(allele, ":")
  short_allele <- lapply(split_allele, function(x, .fields) {
      
    if (length(x) == 0) return(NA) # length(x) is 0 when input is ""
    if (is.na(x[1])) return(NA) 
    
    
    n_fields <- min(.fields, length(x))
    if (n_fields == 1) {
      joined <- paste(x[1], sep = ":")
    }
    if (n_fields == 2) {
      joined <- paste(x[1], x[2], sep = ":")
    }
    if (n_fields == 3) {
      joined <- paste(x[1], x[2], x[3], sep = ":")
    }
    if (n_fields == 4) {
      joined <- paste(x[1], x[2], x[3], x[4], sep = ":")
    }
    return(joined)
  }, .fields = fields)
  
  short_allele <- unlist(short_allele)

  return(short_allele)

}


#' Reduce KIR string to first field
#'
#' Reduce a long KIR string (i.e. "0010101/0020102+0020103/0020104|0030105/0030106/0030107") to keep only
#' the first field, resulting in "001/002+002|003". Keeps all non-numeric entries intact (i.e. "NEG" or "POS")
#'
#' @param kir_string KIR string vector with /, + and | characters
#'
#' @return String vector with first KIR field only, removing duplicates
#'
#' @examples
#' dat <- data.frame(kirstring = c("0010101/0020102+0020103/0020104|0030105/0030106/0030107", "001/002", "NEG", "0010203+0020304|NEG", "POS", "", NA), stringsAsFactors = FALSE)
#' KIR_first_field(dat$kirstring)
KIR_first_field <- function(kir_string) {

  kir_string <- as.character(kir_string)
  
  output <- unlist(lapply(kir_string, function(.kir_string) {
    if (is.na(.kir_string)) {return(NA_character_)}
    if (.kir_string == "") {return(NA_character_)}
    
    # Finds two groups, first group is three digits, followed second group which is any number of
    # digits, and removes the second group
    st <- gsub("(\\d{3})(\\d*)", "\\1", .kir_string)
    
    # Starting from the lowest level, look within each group contained by a + sign for duplicates and remove these, then
    # do the same for the groups contained by a | sign.
    split_or <- strsplit(st, "\\|")
    for (i in 1:length(split_or)) {
      split_plus <- strsplit(split_or[[i]], "\\+")
      for (j in 1:length(split_plus)) {
        split_dash <- strsplit(split_plus[[j]], "/")
        for (k in 1:length(split_dash)) {
          split_dash[[k]] <- paste0(unique(split_dash[[k]]), collapse = "/")
        }
        
        split_plus[[j]] <- paste0(split_dash, collapse = "+")
      }
      
      split_or[[i]] <- paste0(split_plus, collapse = "|")
    }
    
    return(split_or[[1]])
  }))
  
  # TODO: add output verification for no spaces, same number of pluses as in input etc

  return(output)

}
