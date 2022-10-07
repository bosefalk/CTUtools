
#' Load HLA-C Ligand Classification
#'
#' Loads the HLA-C C1/C2 ligand classification database, and keeps the specified first number
#' of allele fields / characters. Removes all "Unclassified" entries in the source data.
#'
#' @param fields Number of fields in the allele to keep (see \code{\link{shorten_allele}}).
#' Default is to keep the first two fields (i.e. 02:05)
#'
#' @return A data.frame with columns "Allele" and "Class"
#' @seealso \code{\link{HLA_C_class_data}}
#'
#' @examples
#' # Default with two fields
#' head(HLA_C_class_load())
#' # With three fields
#' head(HLA_C_class_load(fields = 3))
HLA_C_class_load <- function(fields = 2) {

  #data(HLA_C_class_data)

  HLA_C_class_data$Allele = substring(HLA_C_class_data$Allele,first=3) # Remove C*
  HLA_C_class_data$Allele = sapply(HLA_C_class_data$Allele, FUN = shorten_allele, fields = fields) # from R/munging.R
  HLA_C_class_data = as.data.frame(table(HLA_C_class_data$Allele,HLA_C_class_data$`Predicted Ligand`))
  HLA_C_class_data = HLA_C_class_data[HLA_C_class_data[,3]>=1,]
  HLA_C_class_data = HLA_C_class_data[,-3]
  names(HLA_C_class_data) = c("Allele","Class")
  HLA_C_class_data = subset(HLA_C_class_data, Class != "Unclassified") #alle unclassified raus
  HLA_C_class_data$Class <- as.character(HLA_C_class_data$Class)

  return(HLA_C_class_data)
}


#' Load HLA-B Ligand Classification
#'
#' Loads the HLA-B ligand classification database (Bw6, Bw4 - 80T, Bw4 - 80I classes), and keeps the specified first number
#' of allele fields / characters. Removes all "Unclassified" entries in the source data.
#'
#' @param fields Number of fields in the allele to keep (see \code{\link{shorten_allele}}).
#' Default is to keep the first two fields (i.e. 02:05)
#'
#' @return A data.frame with columns "Allele" and "Class"
#' @seealso \code{\link{HLA_B_class_data}}
#'
#' @examples
#' # Default with two fields
#' head(HLA_B_class_load())
#' # With three fields
#' head(HLA_B_class_load(fields = 3))
HLA_B_class_load <- function(fields = 2) {

  #data(HLA_B_class_data)

  HLA_B_class_data$Allele = substring(HLA_B_class_data$Allele,first=3) # Remove B*
  HLA_B_class_data$Allele = sapply(HLA_B_class_data$Allele, FUN = shorten_allele, fields = fields) # from R/munging.R
  HLA_B_class_data = as.data.frame(table(HLA_B_class_data$Allele,HLA_B_class_data$`Predicted Ligand`))
  HLA_B_class_data = HLA_B_class_data[HLA_B_class_data[,3]>=1,]
  HLA_B_class_data = HLA_B_class_data[,-3]
  names(HLA_B_class_data) = c("Allele","Class")
  HLA_B_class_data = subset(HLA_B_class_data, Class != "Unclassified") #alle unclassified raus
  HLA_B_class_data$Class <- as.character(HLA_B_class_data$Class)

  return(HLA_B_class_data)
}




#' Classify one HLA allele
#'
#' Classifies an allele according to a reference dataframe. If an NMDP code is present, class is given if all possible
#' alleles are the same class, otherwise returns NA
#'
#' @param allele An HLA allele string vector, including ":" characers but not the initial "C*"
#' @param HLA_x_class HLA classification reference data.frame (such as constructed by for example \code{HLA_C_class_load()}), must
#' contain columns Allele (where the number of fields is the same as in \code{allele}) and Class.
#'
#' @return String vector with classification group from reference document, NA if unknown
#'
#' @seealso \code{\link{HLA_C_classification}}, \code{\link{HLA_B_classification}},
#' \code{\link{HLA_C_class_load}}, \code{\link{HLA_B_class_load}}, \code{\link{NMDP}}
#'
#' @examples
#' 
#' dat <- data.frame(HLA_C = c("01:02", "01:AWFCH"), HLA_B = c("07:02", NA),stringsAsFactors = FALSE)
#' HLA_C_class <- HLA_C_class_load()
#' HLA_Classification(dat$HLA_C, HLA_C_class)
#' HLA_B_class <- HLA_B_class_load()
#' HLA_Classification(dat$HLA_B, HLA_B_class)
HLA_Classification = function(allele,HLA_x_class){

  allele <- as.character(allele)
  out <- unlist(lapply(allele, function(xx) {
  #data(NMDP")

  xclass = HLA_x_class$Class[HLA_x_class$Allele==xx] #Wenn kein Analysenfehler und kein NMDP Code
  if (length(xclass)==0) xclass <- NA_character_

  if (xx == "" | is.na(xx)) xclass = NA_character_ # Wenn Analysenfehler

  if (grepl("[a-zA-Z]+",xx)){ #Wenn mit NMDP-Code
    # Classifies all possible alleles that can be constructed with the NMDP code, if they are all the same
    # non-NA value use this classification, otherwise set xclass to NA
    xx = unlist(strsplit(xx, split = ":", fixed = T))
    nmdp_allele = NMDP$Allele[NMDP$nmdp_codes==xx[2]] # TODO: currently only looks at second field
    nmdp_allele = unlist(strsplit(nmdp_allele, split = "/", fixed = T))
    nmdp_allele = gsub("^\\d*:", "", nmdp_allele)
    nmdp_allele = sapply(nmdp_allele, function(x){paste(xx[1],x,sep=":")})
    xclass_split = sapply(nmdp_allele, function(x){
      y = HLA_x_class$Class[HLA_x_class$Allele==x]
      if (length(y) == 0) y = NA
      return(y)
    })

    if (all(is.na(xclass_split))) xclass <- NA_character_
    if (length(unique(xclass_split[!is.na(xclass_split)])) == 1) xclass <- unique(xclass_split[!is.na(xclass_split)])
    else {xclass <- NA_character_}

  }
  return(xclass)
  }))
  
  return(out)
}

#' Classify HLA-C alleles into C1 & C2 groups
#'
#' Takes two HLA-C alleles, classifies them individually into C1 or C2 using \code{\link{HLA_Classification}} and \code{\link{HLA_C_class_data}},
#' and returns the joint C1/C1, C1/C2 or C2/C2 class. If either allel classification is unknown or missing, returns NA.
#'
#' @param allel_c1 Allele column, as string
#' @param allel_c2 Allele column, as string
#' @param fields Number of fields from allele string reference document to use, should be same as number of fields in \code{allel_c1}.
#' For example, if allel_c1 = "02:03", fields should be 2 (the default value). This does not actually apply
#' any string manipulation to \code{allele_c1} & \code{allele_c2}, which needs to be done before passing to this function, 
#' using \code{\link{shorten_allel}}
#' @param HLA_C_class deprecated, used to pre-load classification data before this function was vectorized
#'
#' @return String vector with "C1/C1", "C1/C2", "C2/C2", or \code{NA_character_}
#'
#' @details There are only three outcome classes, C1/C1, C1/C2 and C2/C2. The location of each C-class before or after the "/" does not
#' match the input \code{allele_c1} and \code{_c2} - in other words if \code{allele_c1 = "C2"} and \code{allele_c2 = "C1"} the result is still C1/C2
#' The input alleles should be shortened to match the number of fields in the reference data.
#'
#' @seealso \code{\link{HLA_Classification}}, \code{\link{HLA_C_class_data}}
#'
#' @examples
#' HLA_C_classification("01:02", "01:AWFCH")
#' 
#' dat <- data.frame(C1 = c("01:02", "02:03", "04:10"), C2 = c("01:AWFCH", "07:59", "05:50"))
#' HLA_C_classification(dat$C1, dat$C2)
#' 
#' # If either allele cannot be mapped to the reference data, NA is returned
#' dat_odd <- data.frame(C1 = c("01:02", "01:02", "01:02"), C2 = c(NA, "", "07:02:01"))
#' HLA_C_classification(dat_odd$C1, dat_odd$C2)
HLA_C_classification = function(allele_c1, allele_c2, fields = 2, HLA_C_class = NULL){

  # TODO: Add verification input allele string are same length and structure as nchar
  if (!is.null(HLA_C_class)) {
    warning("HLA_C_class argument is deprecated, classification data is always loaded internally")
  }

  HLA_C_class <- HLA_C_class_load(fields = fields)
  allele_c1 <- as.character(allele_c1)
  allele_c2 <- as.character(allele_c2)
  
  c1 <- unlist(lapply(allele_c1, HLA_Classification, HLA_C_class))
  c2 <- unlist(lapply(allele_c2, HLA_Classification, HLA_C_class))
  tmp_df <- data.frame(c1_class = c1, c2_class = c2, stringsAsFactors = FALSE)
  # Logic for determining joint group, if either allele was not possible to classify returns NA
  out <- dplyr::case_when(
    is.na(tmp_df$c1_class) | is.na(tmp_df$c2_class) ~ NA_character_,
    tmp_df$c1_class == "C2" & tmp_df$c2_class == "C2" ~ "C2/C2",
    tmp_df$c1_class == "C2" & tmp_df$c2_class == "C1" |
      tmp_df$c1_class == "C1" & tmp_df$c2_class == "C2" ~ "C1/C2",
    tmp_df$c1_class == "C1" & tmp_df$c2_class == "C1" ~ "C1/C1",
    TRUE ~ NA_character_
  )

  return(out)
}


#' Classify overall group of HLA-B alleles
#'
#' Takes two HLA-B alleles, classifies them individually into Bw6, Bw4 - 80I or Bw4 - 80T using \code{\link{HLA_Classification}} and
#' \code{\link{HLA_B_class_data}}, and returns the overall group. If either allele is 80I that's the overall group, if not then if either is 80T that's the group,
#' otherwise Bw6.If either allel classification is unknown, returns NA.
#'
#' @param allel_b1 Allele column, as string
#' @param allel_b2 Allele column, as string
#' @param fields Number of fields from allele string reference document to use, should be same as number of fields in \code{allel_b1}.
#' For example, if allel_b1 = "07:02", fields should be 2 (the default value). This does not actually apply
#' any string manipulation to \code{allele_b1} & \code{allele_b2}, which needs to be done before passing to this function, 
#' using \code{\link{shorten_allel}}
#' @param HLA_C_class deprecated, used to pre-load classification data before this function was vectorized
#'
#' @return String vector with "Bw6", "Bw4 - 80T", "Bw4 - 80I", or \code{NA_character_}, as a string
#'
#' @seealso \code{\link{HLA_Classification}}, \code{\link{HLA_B_class_data}}
#'
#' @examples
#' HLA_B_classification("07:02", "08:02")
#' 
#' dat <- data.frame(B1 = c("07:02", "07:02", "07:02"), B2 = c("07:36", "39:15", "08:02"))
#' HLA_B_classification(dat$B1, dat$B2)
#' 
#' # If either allele cannot be mapped to the reference data, NA is returned
#' dat_odd <- data.frame(B1 = c("07:02", "07:02", "07:02"), B2 = c(NA, "", "08:47"))
#' HLA_B_classification(dat_odd$B1, dat_odd$B2)
HLA_B_classification = function(allele_b1, allele_b2, fields = 2, HLA_B_class = NULL){

  # TODO: Add verification input allele string are same length and structure as nchar
  if (!is.null(HLA_B_class)) {
    warning("HLA_B_class argument is deprecated, classification data is always loaded internally")
  }
  
  HLA_B_class <- HLA_B_class_load(fields = fields)
  allele_b1 <- as.character(allele_b1)
  allele_b2 <- as.character(allele_b2)
  
  b1 <- unlist(lapply(allele_b1, HLA_Classification, HLA_B_class))
  b2 <- unlist(lapply(allele_b2, HLA_Classification, HLA_B_class))
  tmp_df <- data.frame(b1_class = b1, b2_class = b2, stringsAsFactors = FALSE)

  # If either is 80I that's the overall group, if not then if either is 80T that's the group, otherwise Bw6.
  out <- dplyr::case_when(
    is.na(tmp_df$b1_class) | is.na(tmp_df$b2_class) ~ NA_character_,
    tmp_df$b1_class == "Bw4 - 80I" | tmp_df$b2_class == "Bw4 - 80I" ~ "Bw4 - 80I",
    tmp_df$b1_class == "Bw4 - 80T" | tmp_df$b2_class == "Bw4 - 80T" ~ "Bw4 - 80T",
    tmp_df$b1_class == "Bw6" & tmp_df$b2_class == "Bw6" ~ "Bw6",
    TRUE ~ NA_character_
  )
  
  return(out)
  
}


#' Assign HLA alleles to Supertype
#' 
#' Assigns a HLA allele to its Supertype (A01, A02 A24 etc). Allele which are classified as Unassigned in the paper are "Unassigned", alleles which do not
#' appear in the paper at all are given "Unknown".
#'
#' @param allele HLA allele string vector as character, two fields only
#' @param HLA "A" for HLA-A, "B" for HLA-B
#'
#' @return string vector with Supertype assignment for each row
#' 
#' @source \url{http://www.biomedcentral.com/1471-2172/9/1}
#' 
#' @examples
#' dat <- data.frame(A_allele = c("01:01", "01:99", NA, "01:13"), 
#'        B_allele = c("07:02", "02:99", NA, "07:10"), stringsAsFactors = FALSE)
#' HLA_Supertype(dat$A_allele, HLA = "A")
#' HLA_Supertype(dat$B_allele, HLA = "B")
HLA_Supertype <- function(allele, HLA) {
  if (!HLA %in% c("A", "B")) {
    stop("HLA input parameter must be either A or B")
  }
  .HLA <- HLA
  super <- Supertype_HLA_lookup %>% filter(HLA == .HLA) %>% select(-HLA)
  
  .allele <- data.frame(allele = as.character(allele), stringsAsFactors = FALSE)
  with_super <- left_join(.allele, super, by = c("allele" = "Allele"))
  
  with_super <- with_super %>% mutate(Supertype = ifelse(is.na(Supertype) & (allele != "" & allele != " " & !is.na(allele)), "Unknown", Supertype))
  
  out <- with_super$Supertype
  return(out)
  
  
  }

