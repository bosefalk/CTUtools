% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HLA_classification.R
\name{HLA_Supertype}
\alias{HLA_Supertype}
\title{Assign HLA alleles to Supertype}
\source{
\url{http://www.biomedcentral.com/1471-2172/9/1}
}
\usage{
HLA_Supertype(allele, HLA)
}
\arguments{
\item{allele}{HLA allele string vector as character, two fields only}

\item{HLA}{"A" for HLA-A, "B" for HLA-B}
}
\value{
string vector with Supertype assignment for each row
}
\description{
Assigns a HLA allele to its Supertype (A01, A02 A24 etc). Allele which are classified as Unassigned in the paper are "Unassigned", alleles which do not
appear in the paper at all are given "Unknown".
}
\examples{
dat <- data.frame(A_allele = c("01:01", "01:99", NA, "01:13"), 
       B_allele = c("07:02", "02:99", NA, "07:10"), stringsAsFactors = FALSE)
HLA_Supertype(dat$A_allele, HLA = "A")
HLA_Supertype(dat$B_allele, HLA = "B")
}
