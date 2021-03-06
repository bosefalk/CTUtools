% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kir.R
\name{KIR3DL1_HLA_B_inhibiting}
\alias{KIR3DL1_HLA_B_inhibiting}
\title{Classifies KIR3DL1 and HLA-B into strong/weak/no inhibiting}
\usage{
KIR3DL1_HLA_B_inhibiting(KIR3DL1_assignment, HLA_B_overall, levels = 2)
}
\arguments{
\item{HLA_B_overall}{String vector with HLA-B Bw classification in format "Bw4 - 80T", from \code{\link{HLA_B_classification}}}

\item{levels}{Can be set to 2, 3 or 4, default is 2}

\item{KIRD3L1_assignment}{String vector, KIR3DL1 assignment in fromat "KIR3DL1-H", from \code{\link{KIR3DL1_3DS1_assignment}}}
}
\value{
If levels = 2, string with "Strong inhibiting", "Weak inhibiting/noninhibiting" or "unknown". If levels = 3, 
split into "Weak inhibiting" and "non-inhibiting". If levels = 4, split "noninhibiting" into "Missing ligand" and "Educated, uninhibiting"
}
\description{
Classifies KIR3DL1 and HLA-B into strong/weak/no inhibiting
}
\examples{
dat <- data.frame(assignment_KIR3DL1 = c("KIR3DL1-L", "KIR3DL1-N", "unknown", "unknown"), HLA_B_group = c("Bw4 - 80T", "Bw4 - 80I", NA, "Bw4 - 80I"), stringsAsFactors = FALSE)
KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group)
KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group, levels = 3)
KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group, levels = 4)
}
