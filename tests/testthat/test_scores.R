context("Boelen scores")

test_that("score_boelen_inhib examples work", {
  dat <- data.frame(kir_2DL1 = "001",
                    kir_2DL2 = "001",
                    kir_2DL3 = c("NEG", "001", "001"),
                    kir_3DL1 = c("NEG", "NEG", "001"),
                    C_class = c("C1/C2", "C2/C2", "C2/C2"),
                    B_class = c("Bw6/Bw6", "Bw6/Bw6", "Bw6/Bw4-80T"),
                    B1 = "07:02",
                    B2 = c("07:02", "46:01", "46:01"),
                    stringsAsFactors = FALSE)
  
  expect_equal(score_boelen_inhib(dat), c(2, 2.75, 2.75))
  expect_equal(score_boelen_inhib(dat, score_or_count = "count"), c(2, 3, 3))
  expect_equal(score_boelen_inhib(dat, include_3DL1 = TRUE), c(2, 2.75, 3.75))
  expect_equal(score_boelen_inhib(dat, separate_2DL2_2DL3 = FALSE), c(2, 2, 2))
  
  dat <- data.frame(kir_2DL1 = "001",
                    kir_2DL2 = "001",
                    kir_2DL3 = c("NEG", "001", "001"),
                    kir_3DL1 = c("NEG", "NEG", "001"),
                    C_class = c("C1/C2", "C2/C2", "C2/C2"),
                    B_class = c("Bw6/Bw6", "Bw6/Bw6", "Bw6/Bw4-80T"),
                    B1 = "07:02",
                    B2 = c("07:02", "46:01", "46:01"),
                    stringsAsFactors = TRUE)
  
  expect_equal(score_boelen_inhib(dat), c(2, 2.75, 2.75))
  expect_equal(score_boelen_inhib(dat, score_or_count = "count"), c(2, 3, 3))
  expect_equal(score_boelen_inhib(dat, include_3DL1 = TRUE), c(2, 2.75, 3.75))
  expect_equal(score_boelen_inhib(dat, separate_2DL2_2DL3 = FALSE), c(2, 2, 2))
  
  
})

test_that("score_boelen_inihb can handle missing and odd values", {
  dat <- data.frame(kir_2DL1 = c(NA, ""),
                    kir_2DL2 = "001",
                    kir_2DL3 = "NEG",
                    kir_3DL1 = "NEG",
                    C_class = "C1/C2",
                    B_class = "Bw6/Bw6",
                    B1 = "07:02",
                    B2 = "07:02",
                    stringsAsFactors = FALSE)
  expect_equal(score_boelen_inhib(dat), c(NA_real_, NA_real_))
  
  dat <- data.frame(kir_2DL1 = "001",
                    kir_2DL2 = c(NA, ""),
                    kir_2DL3 = "NEG",
                    kir_3DL1 = "NEG",
                    C_class = "C1/C2",
                    B_class = "Bw6/Bw6",
                    B1 = "07:02",
                    B2 = "07:02",
                    stringsAsFactors = FALSE)
  expect_equal(score_boelen_inhib(dat), c(NA_real_, NA_real_))
  
  dat <- data.frame(kir_2DL1 = "001",
                    kir_2DL2 = "001",
                    kir_2DL3 = c(NA, ""),
                    kir_3DL1 = "NEG",
                    C_class = "C1/C2",
                    B_class = "Bw6/Bw6",
                    B1 = "07:02",
                    B2 = "07:02",
                    stringsAsFactors = FALSE)
  expect_equal(score_boelen_inhib(dat), c(NA_real_, NA_real_))
  
  dat <- data.frame(kir_2DL1 = "001",
                    kir_2DL2 = "001",
                    kir_2DL3 = "NEG",
                    kir_3DL1 = c(NA, ""),
                    C_class = "C1/C2",
                    B_class = "Bw6/Bw6",
                    B1 = "07:02",
                    B2 = "07:02",
                    stringsAsFactors = FALSE)
  expect_equal(score_boelen_inhib(dat), c(2,2)) # Not using 3DL1 so shouldnt matter that it is missing
  expect_equal(score_boelen_inhib(dat, include_3DL1 = TRUE), c(NA_real_, NA_real_))
  
  dat <- data.frame(kir_2DL1 = "001",
                    kir_2DL2 = "001",
                    kir_2DL3 = "NEG",
                    kir_3DL1 = "NEG",
                    C_class = c(NA, ""),
                    B_class = "Bw6/Bw6",
                    B1 = "07:02",
                    B2 = "07:02",
                    stringsAsFactors = FALSE)
  expect_equal(score_boelen_inhib(dat), c(NA_real_, NA_real_))
  
  dat <- data.frame(kir_2DL1 = "001",
                    kir_2DL2 = "001",
                    kir_2DL3 = "NEG",
                    kir_3DL1 = "NEG",
                    C_class = "C1/C2",
                    B_class = c(NA, ""),
                    B1 = "07:02",
                    B2 = "07:02",
                    stringsAsFactors = FALSE)
  expect_equal(score_boelen_inhib(dat), c(2,2)) # HLA-B ligands only used the 3DL1 are included
  expect_equal(score_boelen_inhib(dat, include_3DL1 = TRUE), c(NA_real_, NA_real_))
  
  dat <- data.frame(kir_2DL1 = "001",
                    kir_2DL2 = "001",
                    kir_2DL3 = "NEG",
                    kir_3DL1 = "NEG",
                    C_class = "C1/C2",
                    B_class = "Bw6/Bw6",
                    B1 = c(NA, ""),
                    B2 = "07:02",
                    stringsAsFactors = FALSE)
  expect_equal(score_boelen_inhib(dat), c(NA_real_, NA_real_))
  
})

