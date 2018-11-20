
#' Load HLA-C Ligand Classification
#'
#' Loads the HLA-C C1/C2 ligand classification database, and keeps the specified first number
#' of allele fields / characters. Removes all "Unclassified" entries in the source data.
#'
#' @param nchar Number of characters in the allele to keep (after removing the inital "C*", but including ":").
#' Default is \code{nchar} = 5, which gives the first two fields (i.e. 02:05)
#'
#' @return A data.frame with columns "Allele" and "Predicted Ligand"
#' @seealso \code{\link{HLA_C_class_data}}
#'
#' @examples
#' # Default with two fields
#' head(HLA_C_class_load())
#' # Only the first field
#' head(HLA_C_class_load(nchar = 2))
HLA_C_class_load <- function(nchar = 5) {

  #data(HLA_C_class_data)

  HLA_C_class_data$Allele = substring(HLA_C_class_data$Allele,first=3) #Ohne C*
  HLA_C_class_data$Allele = substring(HLA_C_class_data$Allele, 1, nchar)
  HLA_C_class_data = as.data.frame(table(HLA_C_class_data$Allele,HLA_C_class_data$`Predicted Ligand`))
  HLA_C_class_data = HLA_C_class_data[HLA_C_class_data[,3]>=1,]
  HLA_C_class_data = HLA_C_class_data[,-3]
  names(HLA_C_class_data) = c("Allele","Class")
  HLA_C_class_data = subset(HLA_C_class_data, Class != "Unclassified") #alle unclassified raus
  # HLA_C_class_data$Class = as.numeric(gsub("C*", "", HLA_C_class_data$`Predicted Ligand`))
  # HLA_C_class_data <- subset(HLA_C_class_data, select = -`Predicted Ligand`)
  # HLA_C_class_data$Class = as.numeric(HLA_C_class_data$Class)
  HLA_C_class_data$Class <- as.character(HLA_C_class_data$Class)

  return(HLA_C_class_data)
}


#' Classify one HLA allele
#'
#' Classifies an allele according to a reference dataframe. If an NMDP code is present, class is given if all possible
#' alleles are the same class, otherwise returns NA
#'
#' @param allel An allele string, including ":" characers but not the initial "C*"
#' @param HLA_x_class HLA classification reference data.frame (such as constructed by for example \code{HLA_C_class_load()})
#'
#' @return String with  classification group from reference document, NA if unknown
#'
#' @seealso \code{\link{HLA_C_class_data}}, \code{\link{NMDP}}
#'
#' @examples
#' HLA_C_class <- HLA_C_class_load()
#' HLA_Classification("01:02", HLA_C_class)
#' HLA_Classification("01:AWFCH", HLA_C_class)
HLA_Classification = function(allele,HLA_x_class){

  #data("data/NMDP.Rdata")

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
#' @param allel_c2 Allele, as string - same length as \code{allel_c1}
#' @param nchar Number of characters from allele string reference document to use, should be same as length as \code{allel_c1}.
#' For example, if allel_c1 = "02:03", nchar should be 5. Default value is 5, for the first two fields.
#'
#' @return One of three C1/C2 groups, as string
#'
#' @seealso \code{\link{HLA_Classification}}, \code{\link{HLA_C_class_data}}, \code{\link{NMDP}}
#'
#' @examples
#' HLA_C_classification("01:02", "01:AWFCH")
HLA_C_classification = function(allele_c1, allele_c2, nchar = 5){

  # TODO: Add verification input allele string are same length and structure as nchar
  HLA_C_class_data <- HLA_C_class_load(nchar = nchar)

  c1 <- HLA_Classification(allele_c1, HLA_C_class_data)
  c2 <- HLA_Classification(allele_c2, HLA_C_class_data)

  # Logic for determining joint group, if either allele was not possible to classify returns NA
  if (any(is.na(c(c1, c2)))) return(NA_character_)
  if (c1 == "C2" & c2 == "C2") return("C2/C2")
  if ((c1 == "C2" & c1 == "C1") | (c1 == "C1" & c1 == "C2")) return("C1/C2")
  if (c1 == "C1" & c2 == "C1") return("C1/C1")
  else {return(NA_character_)}
}

#
# ## HLA-B Klassifizierung
# b1 = HLA_Classification(allel_b1,HLA_B)
# b2 = HLA_Classification(allel_b2,HLA_B)
# HLA_B_class_num = paste(min(b1,b2), max(b1,b2), sep = "/")
# ## Die 11 Proben mit NMDP-Codes bei HLA-B sind ohne Ambiguit?ten klassifizierbar
# #num = HLA_B_class_num
# num = gsub("[0]","",HLA_B_class_num)
# HLA_B_class = as.character(Overall_HLA_B$class[Overall_HLA_B$class_num == num])
# HLA_B_overall = as.character(Overall_HLA_B$overall[Overall_HLA_B$class_num == num])
# if(length(HLA_B_class) == 0){
#   HLA_B_class = "unknown"
#   HLA_B_overall = "unknown"
# }
