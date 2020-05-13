


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
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DL1 & grepl("C2", df$C_class) ~ 1, 
    TRUE ~ 0)  
  
  # 2DL2 is strong functional if together with HLA-C1, B46 or B73, 
  # weak functional with HLA-C2
  df$func_2DL2 <- case_when(
    is.na(df$pres_2DL2) | df$pres_2DL2 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
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
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    is.na(df$B1) | df$B1 == "" ~ NA_real_,
    is.na(df$B2) | df$B2 == "" ~ NA_real_,
    df$pres_2DL3 & (grepl("C1", df$C_class) | 
                                  grepl("^46", df$B1) | grepl("^46", df$B2) | grepl("^73", df$B1) | grepl("^73", df$B2)) ~ score_2dl3, 
    TRUE ~ 0)
  
  # 3DL1 is functional together with Bw4
  if (include_3DL1 == TRUE) {
    df$func_3DL1 <- case_when(
      is.na(df$pres_3DL1) | df$pres_3DL1 == "" ~ NA_real_,
      is.na(df$B_class) | df$B_class == "" | df$B_class == "unknown" ~ NA_real_,
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







#' @noRd
#' 
#' Mapping between commonly used column names from Henning to Bose, i.e. kir2ds1 to kir_2DS1
map_score_df <- function(df) {
  # Henning uses a slightly different naming convention:
  henning_cols <- c("kir2dl1", "kir2dl2", "kir2dl3", "kir3dl1", "kir3dl2", "kir2ds1", "kir2ds2", 
                    "kir2ds4", "kir2ds5", 
                    "hla_c_class.pat", "hla_b_class.pat") 
  try(df <- dplyr::rename(df, kir_2DL1 = kir2dl1), silent = TRUE)
  try(df <- dplyr::rename(df, kir_2DL2 = kir2dl2), silent = TRUE)
  try(df <- dplyr::rename(df, kir_2DL3 = kir2dl3), silent = TRUE)
  try(df <- dplyr::rename(df, kir_3DL1 = kir3dl1), silent = TRUE)
  try(df <- dplyr::rename(df, kir_3DL2 = kir3dl2), silent = TRUE)
  try(df <- dplyr::rename(df, kir_2DS1 = kir2ds1), silent = TRUE)
  try(df <- dplyr::rename(df, kir_2DS2 = kir2ds2), silent = TRUE)
  try(df <- dplyr::rename(df, kir_2DS4 = kir2ds4), silent = TRUE)
  try(df <- dplyr::rename(df, kir_2DS4N = kir2ds4N), silent = TRUE)
  try(df <- dplyr::rename(df, kir_2DS5 = kir2ds5), silent = TRUE)
  try(df <- dplyr::rename(df, kir_3DS1 = kir3ds1), silent = TRUE)
  try(df <- dplyr::rename(df, C_class = hla_c_class.pat), silent = TRUE)
  try(df <- dplyr::rename(df, B_class = hla_b_class.pat), silent = TRUE)
  
  required_cols <- c("kir_2DL1", "kir_2DL2", "kir_2DL3", "kir_3DL1", "kir_3DL2", "kir_2DS1", "kir_2DS2",
                     "kir_2DS4", "kir_2DS5", "kir_3DS1",
                     "C_class", "B_class")
  for (i in 1:length(required_cols)) {
   if (!required_cols[i] %in% names(df)) stop(paste("Column", required_cols[i], "/", henning_cols[i], "not in input dataframe"))
   df[[required_cols[i]]] <- as.character(df[[required_cols[i]]]) # Transform from factor vector
  }
  
  # For HLA allele columns, can either pass a joined column or two separate:
  if (!"hla_a" %in% names(df)) {
    if ("A1" %in% names(df) & "A2" %in% names(df)) {
      df$hla_a <- paste0(df$A1, "/", df$A2) %>% ifelse(grepl("NA", .), NA, .) %>% ifelse(grepl("^/", .) | grepl("/$", .), NA, .)
    }
    if ("a1" %in% names(df) & "a2" %in% names(df)) {
      df$hla_a <- paste0(df$a1, "/", df$a2) %>% ifelse(grepl("NA", .), NA, .) %>% ifelse(grepl("^/", .) | grepl("/$", .), NA, .)
    }
    if (!"hla_a" %in% names(df)) stop("HLA-A column(s) not in input dataframe")
  }
  if (!"hla_b" %in% names(df)) {
    if ("B1" %in% names(df) & "B2" %in% names(df)) {
      df$hla_b <- paste0(df$B1, "/", df$B2) %>% ifelse(grepl("NA", .), NA, .) %>% ifelse(grepl("^/", .) | grepl("/$", .), NA, .)
    }
    if ("b1" %in% names(df) & "b2" %in% names(df)) {
      df$hla_b <- paste0(df$b1, "/", df$b2) %>% ifelse(grepl("NA", .), NA, .) %>% ifelse(grepl("^/", .) | grepl("/$", .), NA, .)
    }
    if (!"hla_b" %in% names(df)) stop("HLA-B column(s) not in input dataframe")
  }
  if (!"hla_c" %in% names(df)) {
    if ("C1" %in% names(df) & "C2" %in% names(df)) {
      df$hla_c <- paste0(df$C1, "/", df$C2) %>% ifelse(grepl("NA", .), NA, .) %>% ifelse(grepl("^/", .) | grepl("/$", .), NA, .)
    }
    if ("c1" %in% names(df) & "c2" %in% names(df)) {
      df$hla_c <- paste0(df$a1, "/", df$a2) %>% ifelse(grepl("NA", .), NA, .) %>% ifelse(grepl("^/", .) | grepl("/$", .), NA, .)
    }
    if (!"hla_c" %in% names(df)) stop("HLA-C column(s) not in input dataframe")
  }
 
  return(df) 
}

#' Krieger inhibitory KIR scores
#' 
#' Calculates the Krieger scores for inhibitory / activating KIRs - see publication below
#'
#' @param data.frame with the following columns:
#' * KIR allele strings, named "kir_2DL1" or "kir2dl1"
#' * HLA-A, B and C allele string, either as one column named "hla_a", or as two "A1" "A2"
#' * C- and B- ligand class column, named "B_class" & "C_class" (or "hla_c_class.pat")
#' @param count_2DS4N_as_2DS4 Default FALSE, if set to TRUE 2DS4 is considered present if 2DS4 or 2DS4N is present
#'
#' @return data.frame with the following columns, all as numeric:
#' * kirl_score
#' * imkir_score
#' * wkir_score
#' 
#' @export
#' @md
#' @source Krieger et al, Killer Immunoglobulin-Like Receptor-Ligand Interactions Predict Clinical Outcomes following Unrelated Donor Transplantations, 
#' Biol Blood Marrow Transplant. 2019 Oct 30, https://www.ncbi.nlm.nih.gov/pubmed/31676338
#'
#' @examples
#' # Create example dataframe
#' dat <- structure(list(kir_2DL1 = 1L, kir_2DL2 = 1L, kir_2DL3 = "NEG", 
#' kir_3DL1 = "NEG", kir_3DL2 = 1L, kir_2DS1 = 1L, kir_2DS2 = "NEG", 
#' kir_2DS4 = 1L, kir_2DS4N = "NEG", kir_2DS5 = 1L, kir_3DS1 = 1L, 
#' C_class = "C1/C1", B_class = "Bw6/Bw4-80T", A1 = "01:01", 
#' A2 = "01:01", B1 = "01:01", B2 = "01:01", C1 = "01:01", C2 = "01:01"), class = "data.frame", row.names = c(NA, 
#'                                                                                                           -1L))
#' score_krieger(dat)
score_krieger = function(df, count_2DS4N_as_2DS4 = FALSE) {
  
  df <- map_score_df(df)
  
  df <- df %>% mutate(pres_2DL1 = CTUtools::KIR_present(kir_2DL1), 
                      pres_2DL2 = CTUtools::KIR_present(kir_2DL2),
                      pres_2DL3 = CTUtools::KIR_present(kir_2DL3),
                      pres_3DL1 = CTUtools::KIR_present(kir_3DL1),
                      pres_3DL2 = CTUtools::KIR_present(kir_3DL2),
                      pres_2DS1 = CTUtools::KIR_present(kir_2DS1),
                      pres_2DS2 = CTUtools::KIR_present(kir_2DS2),
                      pres_2DS4 = CTUtools::KIR_present(kir_2DS4),
                      pres_2DS5 = CTUtools::KIR_present(kir_2DS5),
                      pres_3DS1 = CTUtools::KIR_present(kir_3DS1))
  
  # Count 2DS4 as present if 2DS4N is present
  if (count_2DS4N_as_2DS4 == TRUE) {
    df$pres_2DS4N <- CTUtools::KIR_present(df$kir_2DS4N)
    df$pres_2DS4 <- as.logical(max(df$pres_2DS4, df$pres_2DS4N, na.rm = TRUE))
  }
  
  # For the inhibitory kirs, if the ligand is present the _ligand score is set to 1, and the _missing_ligand score is 0 - if the
  # ligand is missing its the other way around.
  # Activating kirs only record 1 when ligand is present, no missing ligand value
  # These three sets (inhibitory ligand, inhibitory missing ligand, activating ligand) of scores are added together, and
  # used in different configurations for different score calculations - see below
  
  # 2DL1 ligand is C2
  df$kir2dl1_ligand <- case_when(
    is.na(df$pres_2DL1) | df$pres_2DL1 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DL1 & grepl("C2", df$C_class) ~ 1,
    TRUE ~ 0
  )
  df$kir2dl1_missing_ligand <- case_when(
    is.na(df$pres_2DL1) | df$pres_2DL1 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DL1 & !grepl("C2", df$C_class) ~ 1,
    TRUE ~ 0
  )
  
  # 2DL2 ligand is C1
  df$kir2dl2_ligand <- case_when(
    is.na(df$pres_2DL2) | df$pres_2DL2 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DL2 & grepl("C1", df$C_class) ~ 1,
    TRUE ~ 0
  )
  df$kir2dl2_missing_ligand <- case_when(
    is.na(df$pres_2DL2) | df$pres_2DL2 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DL2 & !grepl("C1", df$C_class) ~ 1,
    TRUE ~ 0
  )
  
  # 2DL3 ligand is C1
  df$kir2dl3_ligand <- case_when(
    is.na(df$pres_2DL3) | df$pres_2DL3 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DL3 & grepl("C1", df$C_class) ~ 1,
    TRUE ~ 0
  )
  df$kir2dl3_missing_ligand <- case_when(
    is.na(df$pres_2DL3) | df$pres_2DL3 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DL3 & !grepl("C1", df$C_class) ~ 1,
    TRUE ~ 0
  )
  
  # 3DL1 ligand is Bw4
  df$kir3dl1_ligand <- case_when(
    is.na(df$pres_3DL1) | df$pres_3DL1 == "" ~ NA_real_,
    is.na(df$B_class) | df$B_class == "" | df$B_class == "unknown" ~ NA_real_,
    df$pres_3DL1 & grepl("Bw4", df$B_class) ~ 1,
    TRUE ~ 0
  )
  df$kir3dl1_missing_ligand <- case_when(
    is.na(df$pres_3DL1) | df$pres_3DL1 == "" ~ NA_real_,
    is.na(df$B_class) | df$B_class == "" | df$B_class == "unknown" ~ NA_real_,
    df$pres_3DL1 & !grepl("Bw4", df$B_class) ~ 1,
    TRUE ~ 0
  )
  
  # 3DL2 ligand is A11/A3
  df$kir3dl2_ligand <- case_when(
    is.na(df$pres_3DL2) | df$pres_3DL2 == "" ~ NA_real_,
    is.na(df$hla_a) | df$hla_a == "" ~ NA_real_,
    df$pres_3DL2 & (grepl("11[:]", df$hla_a) | grepl("03[:]", df$hla_a)) ~ 1,
    TRUE ~ 0
  )
  df$kir3dl2_missing_ligand <- case_when(
    is.na(df$pres_3DL2) | df$pres_3DL2 == "" ~ NA_real_,
    is.na(df$hla_a) | df$hla_a == "" ~ NA_real_,
    df$pres_3DL2 & !(grepl("11[:]", df$hla_a) | grepl("03[:]", df$hla_a)) ~ 1,
    TRUE ~ 0
  )
  
  # 2DS1 ligand is C2
  df$kir2ds1_ligand <- case_when(
    is.na(df$pres_2DS1) | df$pres_2DS1 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DS1 & grepl("C2", df$C_class) ~ 1,
    TRUE ~ 0
  )
  
  # 2DS2 ligand is A11
  df$kir2ds2_ligand <- case_when(
    is.na(df$pres_2DS2) | df$pres_2DS2 == "" ~ NA_real_,
    is.na(df$hla_a) | df$hla_a == "" ~ NA_real_,
    df$pres_2DS2 & (grepl("11[:]", df$hla_a)) ~ 1,
    TRUE ~ 0
  )
  
  # 2DS4 ligand is A11
  df$kir2ds4_ligand <- case_when(
    is.na(df$pres_2DS4) | df$pres_2DS4 == "" ~ NA_real_,
    is.na(df$hla_a) | df$hla_a == "" ~ NA_real_,
    df$pres_2DS4 & (grepl("11[:]", df$hla_a)) ~ 1,
    TRUE ~ 0
  )
  
  # 2DS5 ligand is C2
  df$kir2ds5_ligand <- case_when(
    is.na(df$pres_2DS4) | df$pres_2DS4 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DS5 & grepl("C2", df$C_class) ~ 1,
    TRUE ~ 0
  )
  
  
  # ikir is the inhibitory kir score, mkir is missing ligand inhibitory kir score, akir is activating ligand kir score
  df$ikir = df$kir2dl1_ligand+df$kir2dl2_ligand+df$kir2dl3_ligand+df$kir3dl1_ligand+df$kir3dl2_ligand
  df$mkir = df$kir2dl1_missing_ligand+df$kir2dl2_missing_ligand+df$kir2dl3_missing_ligand+
    df$kir3dl1_missing_ligand+df$kir3dl2_missing_ligand
  df$akir = df$kir2ds1_ligand+df$kir2ds2_ligand+df$kir2ds4_ligand+df$kir2ds5_ligand
  
  # 3 possible scores can be constructed from these, these columns are returned from the function
  df$kirl_score = -df$ikir + df$mkir + df$akir
  df$imkir_score = df$ikir + df$mkir
  df$wkir_score = 0.80 * df$ikir + 0.14 * df$akir + 0.99 * df$mkir
  
    return(df %>% select(kirl_score, imkir_score, wkir_score))
}




#' Rafei inhibitory KIR scores
#' 
#' Calculates the Rafei scores for inhibitory / activating KIRs, adjusted version of Krieger scores - see publication below
#'
#' @param df data.frame with the following columns:
#' * KIR allele strings, named "kir_2DL1" or "kir2dl1"
#' * HLA-A, B and C allele string, either as one column named "hla_a", or as two "A1" "A2"
#' * C- and B- ligand class column, named "B_class" & "C_class" (or "hla_c_class.pat")
#' @param count_2DS4N_as_2DS4 Default FALSE, if set to TRUE 2DS4 is considered present if 2DS4 or 2DS4N is present
#'
#' @return data.frame with the following columns, all as numeric:
#' * rafei_inh_kl_matches_2cat
#' * rafei_act_kl_matches_2cat
#' * rafei_inact_kl_matches_2cat
#' 
#' @export
#' @md
#' @source Rafei et al, Role of killer cell immunoglobulin-like receptor (KIR)-ligand interactions to prevent relapse in patients (pts) receiving matched unrelated stem cell transplant (SCT) for acute myeloid leukemia (AML).
#' Journal of Clinical Oncology 37, https://ascopubs.org/doi/abs/10.1200/JCO.2019.37.15_suppl.7049
#'
#' @examples
#' # Create example dataframe
#' dat <- structure(list(kir_2DL1 = 1L, kir_2DL2 = 1L, kir_2DL3 = "NEG", 
#' kir_3DL1 = "NEG", kir_3DL2 = 1L, kir_2DS1 = 1L, kir_2DS2 = "NEG", 
#' kir_2DS4 = 1L, kir_2DS4N = "NEG", kir_2DS5 = 1L, kir_3DS1 = 1L, 
#' C_class = "C1/C1", B_class = "Bw6/Bw4-80T", A1 = "01:01", 
#' A2 = "01:01", B1 = "01:01", B2 = "01:01", C1 = "01:01", C2 = "01:01"), class = "data.frame", row.names = c(NA, 
#'                                                                                                           -1L))
#' score_rafei(dat)
score_rafei <- function(df, count_2DS4N_as_2DS4 = FALSE) {
  
  df <- map_score_df(df)
  
  df <- df %>% mutate(pres_2DL1 = CTUtools::KIR_present(kir_2DL1), 
                      pres_2DL2 = CTUtools::KIR_present(kir_2DL2),
                      pres_2DL3 = CTUtools::KIR_present(kir_2DL3),
                      pres_3DL1 = CTUtools::KIR_present(kir_3DL1),
                      pres_3DL2 = CTUtools::KIR_present(kir_3DL2),
                      pres_2DS1 = CTUtools::KIR_present(kir_2DS1),
                      pres_2DS2 = CTUtools::KIR_present(kir_2DS2),
                      pres_2DS4 = CTUtools::KIR_present(kir_2DS4),
                      pres_2DS5 = CTUtools::KIR_present(kir_2DS5),
                      pres_3DS1 = CTUtools::KIR_present(kir_3DS1))
  
  # Count 2DS4 as present if 2DS4N is present
  if (count_2DS4N_as_2DS4 == TRUE) {
    df$pres_2DS4N <- CTUtools::KIR_present(df$kir_2DS4N)
    df$pres_2DS4 <- as.logical(max(df$pres_2DS4, df$pres_2DS4N, na.rm = TRUE))
  }
  
  # 2DL1 ligand is same as Krieger score, C2
  df$kir2dl1_ligand <- case_when(
    is.na(df$pres_2DL1) | df$pres_2DL1 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    df$pres_2DL1 & grepl("C2", df$C_class) ~ 1,
    TRUE ~ 0
  )

  # 2DL2 ligand v2 is C1 or B46:01 or B73:01 are present
  df$kir2dl2_ligand_v2 <- case_when(
    is.na(df$pres_2DL2) | df$pres_2DL2 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    is.na(df$hla_b) | df$hla_b == "" ~ NA_real_,
    df$pres_2DL2 & grepl("C1", df$C_class) ~ 1,
    df$pres_2DL2 & (grepl("46[:]01",df$hla_b) | grepl("73[:]01",df$hla_b)) ~ 1,
    TRUE ~ 0
  )
  
  # 2DL3 ligand is 2 if both if both C1 and B46:01 or B73:01 are present, 1 if only either C1 or the HLA-B alleles are present and 0 otherwise
  df$kir2dl3_ligand_v2 <- case_when(
    is.na(df$pres_2DL3) | df$pres_2DL3 == "" ~ NA_real_,
    is.na(df$C_class) | df$C_class == "" | df$C_class == "unknown" ~ NA_real_,
    is.na(df$hla_b) | df$hla_b == "" ~ NA_real_,
    df$pres_2DL3 & grepl("C1", df$C_class) & (grepl("46[:]01",df$hla_b) | grepl("73[:]01",df$hla_b)) ~ 2,
    df$pres_2DL3 & grepl("C1", df$C_class) ~ 1,
    df$pres_2DL3 & (grepl("46[:]01",df$hla_b) | grepl("73[:]01",df$hla_b)) ~ 1,
    TRUE ~ 0
  )
  
  # 3DL1 ligand is same as Krieger, Bw4
  df$kir3dl1_ligand <- case_when(
    is.na(df$pres_3DL1) | df$pres_3DL1 == "" ~ NA_real_,
    is.na(df$B_class) | df$B_class == "" | df$B_class == "unknown" ~ NA_real_,
    df$pres_3DL1 & grepl("Bw4", df$B_class) ~ 1,
    TRUE ~ 0
  )
  
  # 3DL2 ligand is same as Krieger, A11/A3
  df$kir3dl2_ligand <- case_when(
    is.na(df$pres_3DL2) | df$pres_3DL2 == "" ~ NA_real_,
    is.na(df$hla_a) | df$hla_a == "" ~ NA_real_,
    df$pres_3DL2 & (grepl("11[:]", df$hla_a) | grepl("03[:]", df$hla_a)) ~ 1,
    TRUE ~ 0
  )
  
  # 2DS1 ligand v2 is no longer just C2, instead only specific HLA-C alleles
  df$kir2ds1_ligand_v2 <- case_when(
    is.na(df$pres_2DS1) | df$pres_2DS1 == "" ~ NA_real_,
    is.na(df$hla_c) ~ NA_real_,
    df$pres_2DS1 & (grepl("02[:]02",df$hla_c) | 
                      grepl("04[:]01",df$hla_c) |
                      grepl("05[:]01",df$hla_c) |
                      grepl("06[:]02",df$hla_c) |
                      grepl("17[:]0",df$hla_c) |
                      grepl("18[:]02",df$hla_c)) ~ 1,
    TRUE ~ 0
  )
  
  # 2DS2 ligand v2 is now specifically A11:01 instead of any A11
  df$kir2ds2_ligand_v2 <- case_when(
    is.na(df$pres_2DS2) | df$pres_2DS2 == "" ~ NA_real_,
    is.na(df$hla_a) | df$hla_a == "" ~ NA_real_,
    df$pres_2DS2 & (grepl("11[:]01", df$hla_a) | grepl("11[:]01", df$hla_a)) ~ 1,
    TRUE ~ 0
  )
  
  # 2DS4 ligand v2 is now a few specific HLA-C alleles as well as A11:01 and A11:02 instead of just A11
  df$kir2ds4_ligand_v2 <- case_when(
    is.na(df$pres_2DS4) | df$pres_2DS4 == "" ~ NA_real_,
    is.na(df$hla_c) | df$hla_c == "" ~ NA_real_,
    is.na(df$hla_a) | df$hla_a == "" ~ NA_real_,
    df$pres_2DS4 & (grepl("01[:]02",df$hla_c) | 
                      grepl("02[:]02",df$hla_c) | 
                      grepl("05[:]01",df$hla_c) |
                      grepl("14[:]02",df$hla_c) |
                      grepl("16[:]01",df$hla_c) |
                      grepl("11[:]01",df$hla_a) |
                      grepl("11[:]02",df$hla_a)) ~ 1,
    TRUE ~ 0
  )
  
  # 2DS5 is no longer usead, instead:
  # 3DS1 ligant is Bw4-80T or HLA-B 27:05
  df$kir3ds1_ligand <- case_when(
    is.na(df$pres_3DS1) | df$pres_3DS1 == "" ~ NA_real_,
    is.na(df$B_class) | df$B_class == "" | df$B_class == "unknown" ~ NA_real_,
    is.na(df$hla_b) | df$hla_b == "" ~ NA_real_,
    df$pres_3DS1 & grepl("80T", df$B_class) ~ 1,
    df$pres_3DS1 & grepl("27[:]05",df$hla_b) ~ 1,
    TRUE ~ 0
  )
  
 df$ikir_score_v2 = df$kir2dl1_ligand + df$kir2dl2_ligand_v2 + df$kir2dl3_ligand_v2 + 
   df$kir3dl1_ligand + df$kir3dl2_ligand
 df$akir_score_v2 = df$kir2ds1_ligand_v2 + df$kir2ds2_ligand_v2 + df$kir2ds4_ligand_v2 + df$kir3ds1_ligand
 
 
 df$rafei_inh_kl_matches_2cat = case_when(
   is.na(df$ikir_score_v2) | is.na(df$akir_score_v2) ~ NA_character_,
   df$ikir_score_v2 >= 3 ~ ">=3", 
   TRUE ~ "<3")
 df$rafei_act_kl_matches_2cat = case_when(
   is.na(df$ikir_score_v2) | is.na(df$akir_score_v2) ~ NA_character_,
   df$akir_score_v2 >= 1 ~ ">=1",
   TRUE ~ "0"
 )
   
  df$rafei_inact_kl_matches_2cat = case_when(
    is.na(df$ikir_score_v2) | is.na(df$akir_score_v2) ~ NA_character_,
    (df$ikir_score_v2 >= 3 & df$akir_score_v2 == 0) ~ "unfav",
    TRUE ~ "fav"
  )
    
 return(df %>% select(starts_with("rafei")))
 
}

