#' @section Main functions:
#'
#' \code{\link{shorten_allele}}: Reduce a long allele string (i.e. "01:02:01:03g") down to fewer fields ("01:02")
#'
#' \code{\link{HLA_C_classification}}: Classify two HLA-C alleles into a C1/C2 group
#'
#' \code{\link{HLA_B_classification}}: Classify two HLA-B alleles into an overall group Bw6, Bw4-80I or -80T
#'
#' \code{\link{KIR_det_GCN}}: Convert a KIR string (such as "004+003|006+010") into a Gene Copy Number
#'
#'
#' @section Datasets:
#'
#' \code{\link{HLA_C_class_data}}: Reference data for C1 / C2 classification
#'
#' \code{\link{HLA_B_class_data}}: Reference data for Bw6, Bw4-80I and Bw4-80T classification
#'
#' \code{\link{NMDP}}: Reference data for expanding NMDP codes
#'
#' @section Location:
#'
#' S:/Donor/Rlib
#'
#' \code{packrat::set_opts(local.repos = "S://Donor//Rlib")}
#'
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"
