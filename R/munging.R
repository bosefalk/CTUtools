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
