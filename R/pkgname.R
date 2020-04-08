#' @section Main functions:
#'
#' \code{\link{shorten_allele}}: Reduce a long allele string (i.e. "01:02:01:03g") down to fewer fields ("01:02")
#'
#' \code{\link{KIR_first_field}}: Reduce a long KIR string down to only use the first field (i.e. "001+002|003")
#'
#' \code{\link{KIR_present}}: Determine if a KIR is present or absent given a KIR string
#'
#' \code{\link{HLA_C_classification}}: Classify two HLA-C alleles into a C1/C2 group
#'
#' \code{\link{HLA_B_classification}}: Classify two HLA-B alleles into an overall group Bw6, Bw4-80I or -80T
#'
#' \code{\link{KIR_det_GCN}}: Convert a KIR string (such as "004+003|006+010") into a Gene Copy Number
#' 
#' \code{\link{KIR3DL1_HLA_B_inhibiting}}: Classifies KIR3DL1 and HLA-B into strong/weak/no inhibiting
#' 
#' \code{\link{KIR3DL1_3DS1_assignment}}: Classifies KIR 3DL1 & 3DS1 alleles into KIR3DL1-H, -L, -N
#' 
#' @section Score functions:
#' 
#' \code{\link{score_boelen_inhib}}: Boelen inhibitory score calculation
#' 
#' \code{\link{score_krieger}}: Krieger inhibitory and activating KIR score calculation
#' 
#' \code{\link{score_rafei}}: Rafei varation of Kreiger scores
#'
#' @section Datasets:
#'
#' \code{\link{HLA_C_class_data}}: Reference data for C1 / C2 classification
#'
#' \code{\link{HLA_B_class_data}}: Reference data for Bw6, Bw4-80I and Bw4-80T classification
#'
#' \code{\link{NMDP}}: Reference data for expanding NMDP codes
#'
#' \code{\link{ASSIGN_KIR3DL1}}: KIR3DL1/HLA-B Subtypes, from Bodreau et al
#'
#' @section Location:
#'
#' \url{https://github.com/bosefalk/CTUtools}
#'
#' \code{devtools::install_github("bosefalk/CTUtools")}
#'
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"
