#' Reduce allele string to specified number of fields
#'
#' Reduce a long allele string (i.e. "01:02:01:03g") down to fewer fields,
#' default being 2 fields. NMDP codes are preserved, trailing "g" are removed.
#'
#' @param allele String, including ":" but not initial "A*", so for example "01:02:01:03"
#' @param fields Number of fields to keep, default = 2 which would return "01:02"
#'
#' @return String with shorterned fields
#'
#' @examples
#' shorten_allele("01:02:01:03")
#' shorten_allele("03:ADJRE")
#' shorten_allele("01:02:01:03", fields = 3)
shorten_allele <- function(allele, fields = 2) {


  if (!fields %in% c(1,2,3,4)) {
    stop("shortern_allele: fields can only be set to 1,2,3 or 4")
  }



  # Remove trailing g if present
  allele <- gsub("G$", "", allele)

  split_allele <- strsplit(allele, ":")[[1]]
  if (fields == 1) {
    short_allele <- paste0(split_allele[1])
  }

  if (fields == 2) {
    short_allele <- paste0(split_allele[1], ":", split_allele[2])
  }

  #TODO: Check this works for NMDP fields
  if (fields == 3) {
    short_allele <- paste0(split_allele[1], ":", split_allele[2], ":", split_allele[3])
  }

  #TODO: Check this works for NMDP fields
  if (fields == 4) {
    short_allele <- paste0(split_allele[1], ":", split_allele[2], ":", split_allele[3], ":", split_allele[4])
  }

  return(short_allele)

}

#' Reduce KIR string to first field
#'
#' Reduce a long KIR string (i.e. "0010101/0020102+0020103/0020104|0030105/0030106/0030107") to keep only
#' the first field, resulting in "001/002+002|003". Keeps all non-numeric entries intact (i.e. "NEG" or "POS")
#'
#' @param kir_string KIR string with /, + and | characters
#'
#' @return String with first KIR field only, removing duplicates
#'
#' @examples
#' KIR_first_field("0010101/0020102+0020103/0020104|0030105/0030106/0030107")
#' KIR_first_field("001/002")
#' KIR_first_field("0010203+0020304|NEG")
#' KIR_first_field("POS")
KIR_first_field <- function(kir_string) {

  # Finds two groups, first group is three digits, followed second group which is any number of
  # digits, and removes the second group
  st <- gsub("(\\d{3})(\\d*)", "\\1", kir_string)

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

  st <- split_or

  # TODO: add output verification for no spaces, same number of pluses as in input etc

  return(st[[1]])

}
