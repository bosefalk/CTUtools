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
#' @seealso \code{\link{HLA_C_class_load}}
#' @source \url{https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?C}
"HLA_C_class_data"

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

