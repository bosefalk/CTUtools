.onLoad <- function(libname, pkgname) {
  data("HLA_C_class_data", "HLA_B_class_data", "NMDP",
       package=pkgname, envir=parent.env(environment()))
}
