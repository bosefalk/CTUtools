% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{score_krieger}
\alias{score_krieger}
\title{Krieger inhibitory KIR scores}
\source{
Krieger et al, Killer Immunoglobulin-Like Receptor-Ligand Interactions Predict Clinical Outcomes following Unrelated Donor Transplantations,
Biol Blood Marrow Transplant. 2019 Oct 30, https://www.ncbi.nlm.nih.gov/pubmed/31676338
}
\usage{
score_krieger(df, count_2DS4N_as_2DS4 = FALSE)
}
\arguments{
\item{count_2DS4N_as_2DS4}{Default FALSE, if set to TRUE 2DS4 is considered present if 2DS4 or 2DS4N is present}

\item{data.frame}{with the following columns:
\itemize{
\item KIR allele strings, named "kir_2DL1" or "kir2dl1"
\item HLA-A, B and C allele string, either as one column named "hla_a", or as two "A1" "A2"
\item C- and B- ligand class column, named "B_class" & "C_class" (or "hla_c_class.pat")
}}
}
\value{
data.frame with the following columns, all as numeric:
\itemize{
\item kirl_score
\item imkir_score
\item wkir_score
}
}
\description{
Calculates the Krieger scores for inhibitory / activating KIRs - see publication below
}
\examples{
# Create example dataframe
dat <- structure(list(kir_2DL1 = 1L, kir_2DL2 = 1L, kir_2DL3 = "NEG", 
kir_3DL1 = "NEG", kir_3DL2 = 1L, kir_2DS1 = 1L, kir_2DS2 = "NEG", 
kir_2DS4 = 1L, kir_2DS4N = "NEG", kir_2DS5 = 1L, kir_3DS1 = 1L, 
C_class = "C1/C1", B_class = "Bw6/Bw4-80T", A1 = "01:01", 
A2 = "01:01", B1 = "01:01", B2 = "01:01", C1 = "01:01", C2 = "01:01"), class = "data.frame", row.names = c(NA, 
                                                                                                          -1L))
score_krieger(dat)
}
