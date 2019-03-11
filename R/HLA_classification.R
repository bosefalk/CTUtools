
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
#' @param allele An allele string, including ":" characers but not the initial "C*"
#' @param HLA_x_class HLA classification reference data.frame (such as constructed by for example \code{HLA_C_class_load()}), must
#' contain columns Allele (where the number of fields is the same as in \code{allele}) and Class.
#'
#' @return String with classification group from reference document, NA if unknown
#'
#' @seealso \code{\link{HLA_C_classification}}, \code{\link{HLA_B_classification}},
#' \code{\link{HLA_C_class_load}}, \code{\link{HLA_B_class_load}}, \code{\link{NMDP}}
#'
#' @examples
#' HLA_C_class <- HLA_C_class_load()
#' HLA_Classification("01:02", HLA_C_class)
#' HLA_Classification("01:AWFCH", HLA_C_class)
#' HLA_B_class <- HLA_B_class_load()
#' HLA_Classification("07:02", HLA_B_class)
HLA_Classification = function(allele,HLA_x_class){

  #data(NMDP")

  xx = allele

  xclass = HLA_x_class$Class[HLA_x_class$Allele==xx] #Wenn kein Analysenfehler und kein NMDP Code
  if (length(xclass)==0) xclass <- NA_character_

  if (xx == "" | is.na(xx)) xclass = NA_character_ # Wenn Analysenfehler

  if (grepl("[a-zA-Z]+",xx)){ #Wenn mit NMDP-Code
    # Classifies all possible alleles that can be constructed with the NMDP code, if they are all the same
    # non-NA value use this classification, otherwise set xclass to NA
    xx = unlist(strsplit(xx, split = ":", fixed = T))
    nmdp_allele = NMDP$Allele[NMDP$nmdp_codes==xx[2]] # TODO: currently only looks at second field
    nmdp_allele = unlist(strsplit(nmdp_allele, split = "/", fixed = T))
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
}

#' Classify HLA-C alleles into C1 & C2 groups
#'
#' Takes two HLA-C alleles, classifies them individually into C1 or C2 using \code{\link{HLA_Classification}} and \code{\link{HLA_C_class_data}},
#' and returns the joint C1/C1, C1/C2 or C2/C2 class. If either allel classification is unknown, returns NA.
#'
#' @param allel_c1 Allele, as string
#' @param allel_c2 Allele, as string
#' @param fields Number of fields from allele string reference document to use, should be same as number of fields in \code{allel_c1}.
#' For example, if allel_c1 = "02:03", fields should be 2 (the default value). This does not actually apply
#' any string manipulation to \code{allele_c1} & \code{c2}, which needs to be done before passing to this function.
#' @param HLA_C_class optional, dataframe from \code{\link{HLA_C_class_load}}, if not supplied this will be called within the function
#' This option is here to pre-load this reference data, in case where you need to apply this function across a large number of subjects.
#'
#' @return One of "C1/C1", "C1/C2", "C2/C2", or \code{NA_character_}, as a string
#'
#' @details There are only three outcome classes, C1/C1, C1/C2 and C2/C2. The location of each C-class before or after the "/" does not
#' match the input \code{allele_c1} and \code{_c2} - in other words if \code{allele_c1 = "C2"} and \code{allele_c2 = "C1"} the result is still C1/C2
#'
#' @seealso \code{\link{HLA_Classification}}, \code{\link{HLA_C_class_data}}
#'
#' @examples
#' HLA_C_classification("01:02", "01:AWFCH")
#'
#' # For when you need to optimize execution time, pre-load the reference database
#' referencedata <- HLA_C_class_load()
#' HLA_C_classification("01:02", "01:AWFCH", HLA_C_class = referencedata)
HLA_C_classification = function(allele_c1, allele_c2, fields = 2, HLA_C_class = NULL){

  # TODO: Add verification input allele string are same length and structure as nchar
  if (is.null(HLA_C_class)) {
    HLA_C_class <- HLA_C_class_load(fields = fields)
  }

  c1 <- HLA_Classification(allele_c1, HLA_C_class)
  c2 <- HLA_Classification(allele_c2, HLA_C_class)

  # Logic for determining joint group, if either allele was not possible to classify returns NA
  if (any(is.na(c(c1, c2)))) return(NA_character_)
  if (c1 == "C2" & c2 == "C2") return("C2/C2")
  if ((c1 == "C2" & c2 == "C1") | (c1 == "C1" & c2 == "C2")) return("C1/C2")
  if (c1 == "C1" & c2 == "C1") return("C1/C1")
  else {return(NA_character_)}
}


#' Classify overall group of HLA-B alleles
#'
#' Takes two HLA-B alleles, classifies them individually into Bw6, Bw4 - 80I or Bw4 - 80T using \code{\link{HLA_Classification}} and
#' \code{\link{HLA_B_class_data}}, and returns the overall group. If either allele is 80I that's the overall group, if not then if either is 80T that's the group,
#' otherwise Bw6.If either allel classification is unknown, returns NA.
#'
#' @param allele_b1 Allele, as string
#' @param allele_b2 Allele, as string
#' @param fields Number of fields from allele string reference document to use, should be same as number of fields in \code{allel_b1}.
#' For example, if allel_b1 = "44:180", fields should be 2 (the default value). This does not actually apply
#' any string manipulation to \code{allele_b1} & \code{b2}, which needs to be done before passing to this function.
#' @param HLA_B_class optional, dataframe from \code{\link{HLA_B_class_load}}, if not supplied this will be called within the function
#' This option is here to pre-load this reference data, in case where you need to apply this function across a large number of subjects.
#'
#'
#' @return One of "Bw6", "Bw4 - 80T", "Bw4 - 80I", or \code{NA_character_}, as a string
#'
#' @seealso \code{\link{HLA_Classification}}, \code{\link{HLA_B_class_data}}
#'
#' @examples
#' HLA_B_classification("44:180", "07:02")
#' # For when you need to optimize execution time, pre-load the reference database
#' referencedata <- HLA_B_class_load()
#' HLA_B_classification("44:180", "07:02", HLA_B_class = referencedata)
HLA_B_classification = function(allele_b1, allele_b2, fields = 2, HLA_B_class = NULL){

  # TODO: Add verification input allele string are same length and structure as nchar
  if (is.null(HLA_B_class)) {
    HLA_B_class <- HLA_B_class_load(fields = fields)
  }
  
  b1 = HLA_Classification(allele_b1,HLA_B_class)
  b2 = HLA_Classification(allele_b2,HLA_B_class)

  # If either is 80I that's the overall group, if not then if either is 80T that's the group, otherwise Bw6.
  if (any(is.na(c(b1, b2)))) {return(NA_character_)}
  if (any(c(b1, b2) == "Bw4 - 80I")) {return("Bw4 - 80I")}
  if (any(c(b1, b2) == "Bw4 - 80T")) {return("Bw4 - 80T")}
  if (all(c(b1, b2) == "Bw6")) {return("Bw6")}

}
