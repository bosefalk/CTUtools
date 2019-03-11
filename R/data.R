#' HLA-C C1/C2 ligand classification database
#'
#' Classification database of HLA-C Ligands from Alleles, from Immuno Polymorphism Database
#'
#'
#' @format
#' \describe{
#'   \item{Allele}{HLA-C allele, in format "C*01:02:01:04"}
#'   \item{Predicted Ligand}{Predicted Ligand, "C1" or "C2"}
#' }
#'
#' @seealso \code{\link{HLA_C_class_load}}, \code{\link{HLA_C_classification}}
#' @source \url{https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?C}
"HLA_C_class_data"

#' HLA-B Bw4/Bw6-80I/T ligand classification database
#'
#' Classification database of HLA-B Ligands from Alleles, from Immuno Polymorphism Database
#'
#'
#' @format
#' \describe{
#'   \item{Allele}{HLA-B allele, in format "B*01:02:01:04"}
#'   \item{Predicted Ligand}{Predicted Ligand, "Bw6", "Bw4 - 80I" or "Bw4 - 80T"}
#' }
#'
#' @seealso \code{\link{HLA_B_class_load}}, \code{\link{HLA_B_classification}}
#' @source \url{https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?B}
"HLA_B_class_data"


#' NMDP allele codes
#'
#' Lists possible alleles for NMDP codes
#'
#'
#' @format
#' \describe{
#'   \item{nmdp_codes}{NMDP code}
#'   \item{Allele}{Component alleles, separated by "/"}
#' }
#'
#' @source \url{https://hml.nmdp.org/MacUI/}
"NMDP"

#' KIR3DL1/HLA-B Subtypes
#'
#' KIR3DL1/HLA-B Subtypes Govern Acute Myelogenous Leukemia Relapse After Hematopoietic Cell Transplantation, Boudreau et al. 
#'
#' @format
#' \describe{
#'   \item{Allele}{KIR3DL1 allele}
#'   \item{Class}{Allele classification. 0 means unclassified}
#' }
#' @seealso \code{\link{KIR3DL1_3DS1_assignment}}
#' @source Boudreau et al. Journal of Clinical Oncology 2017 35:20, 2268-2278
"ASSIGN_KIR3DL1"
