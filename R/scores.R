


#' Boelen inhibitory score
#' 
#' Calculates the Boelen inhibitory score for a set of KIR 2DL1, 2DL2, 2DL3, HLA-B and -C ligands and HLA-B alleles
#' (optionally also KIR 3DL1). Several options for how the score is calculated - see Details and parameter inputs
#'
#' @param df data frame with the following columns - if any of these columns are not present with these names, the function will stop:
#' * kir_2DL1, kir_2DL2, kir_2DL3, kir_3DL1 all in standard KIR string format ("001+002", "NEG").
#' * Columns C_class and B_class are HLA-C and -B ligands, in format "C1/C1" and "Bw6/Bw4-80I" - see \code{\link{HLA_C_classification}} and \code{\link{HLA_B_classification}}.
#' * HLA-B alleles should be in two columns B1 and B2, standard format "07:02"
#' @param include_3DL1 default FALSE, if KIR 3DL1 score is included in calculation
#' @param separate_2DL2_2DL3 default TRUE, whether to add the scores from 2DL2 and 2DL3 separately to the total, or if
#' the max of the two should be used
#' @param score_or_count default "score", can be set to "count" - see details
#'
#' @details   The score is calculated by:
#' (1 if Functional 2DL1) + (1 if Strong Functional 2DL2 or 0.5 if weak Functional 2DL2) + 
#' (0.75 if Functional 2DL3) + (1 if Functional 3DL1)
#' 
#'   If instead the count is given:
#'   Inhibitory count = 1 for each of Functional 2DL1, Functional 2DL2/L3, Functional 3DL1
#'
#' What constitues a "Functional" KIR differs:
#' * 2DL1 is functional if together with HLA-C2 ligand
#' * 2DL2 is strong functional if together with HLA-C1 ligand, HLA-B alleles B46 or B73, weak functional with HLA-C2 ligand
#' * 2DL3 is functional if together with HLA-C1 ligand, HLA-B alleles B46 or B73
#' * 3DL1 is functional together with Bw4
#' @md
#' @return Vector of scores for each row in the input dataset
#' 
#' @source Boudreau et al. Science Immunology 09 Nov 2018: Vol. 3, Issue 29, eaao2892 DOI: 10.1126/sciimmunol.aao2892
#'
#' @examples
#' dat <- data.frame(kir_2DL1 = "001",
#'                   kir_2DL2 = "001",
#'                   kir_2DL3 = c("NEG", "001", "001"),
#'                   kir_3DL1 = c("NEG", "NEG", "001"),
#'                   C_class = c("C1/C2", "C2/C2", "C2/C2"),
#'                   B_class = c("Bw6/Bw6", "Bw6/Bw6", "Bw6/Bw4-80T"),
#'                   B1 = "07:02",
#'                   B2 = c("07:02", "46:01", "46:01"),
#'                   stringsAsFactors = FALSE)
#'                   
#' score_boelen_inhib(dat)
#' score_boelen_inhib(dat, score_or_count = "count")
#' score_boelen_inhib(dat, include_3DL1 = TRUE)
#' score_boelen_inhib(dat, separate_2DL2_2DL3 = FALSE)
score_boelen_inhib <- function(df, include_3DL1 = FALSE, separate_2DL2_2DL3 = TRUE, 
                         score_or_count = "score") {
 
  # For each subject, need info on kir presence / absence for 2DL1, L2, L3 and optionally 3DL1 - this is found
  # from the raw KIR string in columns kir_2DL1 etc (i.e. "001+003", "NEG")
  # HLA-C class ligand is one of "C1/C1", "C1/C2" or "C2/C2" in C_class column
  # Raw HLA-B allele info is also needed, in columns B1 and B2 (i.e. "02:07", "43:01")
  # When including 3DL1, also needed is the B-ligand classification, in B_class_B1 and B_class_B2 ("Bw4-80T", "Bw6")
  # this function assumes the dataframe has these exact column names, need to do data munging beforehand to make df look like this
  
  required_cols <- c("kir_2DL1", "kir_2DL2", "kir_2DL3",
                     "C_class", "B_class", "B1", "B2")
  if (include_3DL1 == TRUE) {required_cols <- append(required_cols, "kir_3DL1")}
  for (l in required_cols) {
    if (!l %in% names(df)) stop(paste("Column", l, "not in input dataframe"))
    df[[l]] <- as.character(df[[l]]) # Transform from vector
  }
  
  df <- df %>% mutate(pres_2DL1 = CTUtools::KIR_present(kir_2DL1), 
                      pres_2DL2 = CTUtools::KIR_present(kir_2DL2),
                      pres_2DL3 = CTUtools::KIR_present(kir_2DL3),
                      pres_3DL1 = CTUtools::KIR_present(kir_3DL1))
  
  # Split B_Class into components
  df = df %>% mutate(B_class_B1 = sapply(df$B_class, function(x) unlist(strsplit(x, "/"))[1])) %>% 
    mutate(B_class_B2 = sapply(df$B_class, function(x) unlist(strsplit(x, "/"))[2]))
  
  # Score is given by:
  # (1 if Func 2DL1) + (1 if Strong Func 2DL2 or 0.5 if weak Func 2DL2) + (0.75 if Func 2DL3) + (1 if Func 3DL1)
  
  # Count is given by:
  # Inhibitory count = 1 for each of Func 2DL1, Func 2DL2/L3, Func 3DL1
  
  score_2dl2 = ifelse(score_or_count == "score", .5, 1)
  score_2dl3 = ifelse(score_or_count == "score", .75, 1)
  
  # 2DL1 is functional if together with HLA-C2
  df$func_2DL1 <- case_when(
    is.na(df$pres_2DL1) | df$pres_2DL1 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" ~ NA_real_,
    df$pres_2DL1 & grepl("C2", df$C_class) ~ 1, 
    TRUE ~ 0)  
  
  # 2DL2 is strong functional if together with HLA-C1, B46 or B73, 
  # weak functional with HLA-C2
  df$func_2DL2 <- case_when(
    is.na(df$pres_2DL2) | df$pres_2DL2 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" ~ NA_real_,
    is.na(df$B1) | df$B1 == "" ~ NA_real_,
    is.na(df$B2) | df$B2 == "" ~ NA_real_,
    df$pres_2DL2 & (grepl("C1", df$C_class) | 
                                  grepl("^46", df$B1) | grepl("^46", df$B2) | grepl("^73", df$B1) | grepl("^73", df$B2)) ~ 1,
    df$pres_2DL2 & grepl("C2", df$C_class) ~ score_2dl2, # case_when is heirarcical, so will hit C1 condition and reutn 1 before checking this
    TRUE ~ 0
  )
  
  # 2DL3 is functional if together with HLA-C1, B46 or B73
  df$func_2DL3 <- case_when(
    is.na(df$pres_2DL3) | df$pres_2DL3 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" ~ NA_real_,
    is.na(df$B1) | df$B1 == "" ~ NA_real_,
    is.na(df$B2) | df$B2 == "" ~ NA_real_,
    df$pres_2DL3 & (grepl("C1", df$C_class) | 
                                  grepl("^46", df$B1) | grepl("^46", df$B2) | grepl("^73", df$B1) | grepl("^73", df$B2)) ~ score_2dl3, 
    TRUE ~ 0)
  
  # 3DL1 is functional together with Bw4
  if (include_3DL1 == TRUE) {
    df$func_3DL1 <- case_when(
      is.na(df$pres_3DL1) | df$pres_3DL1 == "" ~ NA_real_,
      is.na(df$B_class) | df$B_class == "" ~ NA_real_,
      df$pres_3DL1 & (grepl("Bw4", df$B_class_B1) | grepl("Bw4", df$B_class_B2)) ~ 1, 
      TRUE ~ 0)
  }
  
  
  df$score <- df$func_2DL1
  
  if (separate_2DL2_2DL3 == FALSE) { # If 2DL2 & L3 are counted together just use the max of the two values
    df$score <- df$score + pmax(df$func_2DL2, df$func_2DL3)
  } else {
    df$score <- df$score + df$func_2DL2 + df$func_2DL3
  }
  
  
  if (include_3DL1 == TRUE) {
    df$score <- df$score + df$func_3DL1
  }
  
  
  return(df$score)
}


