#' Determine KIR Gene Copy Number
#'
#' @param kir_string String with KIR info with the first field only (three characters), example "004+003|006+010"
#' @param return_numeric If set to TRUE, returns numeric GCN instead of character, sets "x" to NA and 2|3 to 2.5 etc
#'
#' @return Default a string with "x" if "POS", "0" if NEG, NA if NA or empty string, and otherwise the GCN. If \code{return_numeric = TRUE}
#' instead returns a numeric value, where x results are NA and "2|3" are 2.5
#'
#' @examples
#'
#' str(KIR_det_GCN("004+003|006+010"))
#' str(KIR_det_GCN("004+003|006+010", return_numeric = TRUE))
#' KIR_det_GCN("POS")
#' KIR_det_GCN("POS", return_numeric = TRUE)
KIR_det_GCN = function(kir_string, return_numeric = FALSE){
  require(dplyr)

  # Call this on the KIR string if it's not empty/NA/POS/NEG:
  det_clean_GCN <- function(y) {
    split = unlist(strsplit(y,split = "[|]"))
    clean_gcn = nchar(split) - nchar(gsub("[+]","",split))+1
    clean_gcn = ifelse(min(clean_gcn) == max(clean_gcn),min(clean_gcn),paste(min(clean_gcn),max(clean_gcn),sep="|"))
    return(as.character(clean_gcn))
  }

  # Checks conditions according to top first, then working down and applying GCN determination function on string if none of
  # the other conditions are fulfilled
  gcn <- dplyr::case_when(
    grepl("POS", kir_string) ~ "x",
    grepl("NEG", kir_string) ~ "0",
    is.na(kir_string) | grepl("^$", kir_string) ~ NA_character_,
    TRUE ~ det_clean_GCN(kir_string)
  )

  if (return_numeric == TRUE) {
    # When running src/KIR_investigation.Rmd from 18-01 project error is when hitting gcn==1|2 statement no TRUE/FALSE value
    try({
      if(is.na(gcn)) {gcn <- NA_integer_}
      if(gcn == "1|2") {gcn <- 1.5}
      if(gcn == "2|3") {gcn <- 2.5}
      if(gcn == "3|4") {gcn <- 3.5}
      gcn <- suppressWarnings(as.numeric(gcn))}, silent = TRUE)
    return(gcn)
  } else {return(gcn)}
}
