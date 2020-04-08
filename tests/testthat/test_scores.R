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

context("score mapping dataframe function")

mapdf_dat <- data.frame(kir2dl1 = "001", 
                        kir2dl2 = "001",
                        kir2dl3 = "001",
                        kir3dl1 = "001",
                        kir3dl2 = "001",
                        kir2ds1 = "001",
                        kir2ds2 = "001",
                        kir2ds4 = "001",
                        kir2ds5 = "001",
                        kir3ds1 = "001",
                        hla_c_class.pat = "C1/C1",
                        hla_b_class.pat = "Bw6/Bw6",
                        hla_a = "01:01/01:01", 
                        hla_b = "01:01/01:01", 
                        hla_c = "01:01/01:01")
required_columns <- c("kir_2DL1", "kir_2DL2", "kir_2DL3", "kir_3DL1", "kir_3DL2", 
                      "kir_2DS1", "kir_2DS2", "kir_2DS4", "kir_2DS5", "C_class", "B_class", 
                      "hla_a", "hla_b", "hla_c")
test_that("score mapping function works for transforming Hennings column names", {
  expect_true(all(required_columns %in% names(map_score_df(mapdf_dat))))
  dat2 <- mapdf_dat %>% select(-kir2ds4)
  expect_error(mapdf_dat(dat2))
})

test_that("score mapping function correctly deals with 2DS4N", {
  dat3 <- mapdf_dat %>% mutate(kir2ds4N = "NEG")
  expect_true("kir_2DS4N" %in% names(map_score_df(dat3)))
})

test_that("score mapping HLA columns different combinations", {
  dat <- mapdf_dat %>% select(-hla_a, -hla_b, -hla_c) %>% 
    mutate(A1 = 1, A2 = 1, B1 = 1, B2 = 2, C1 = 1, C2 = 1)
  out <- map_score_df(dat)
  expect_true(all(c("hla_a", "hla_b", "hla_c") %in% names(out)))
  expect_true(!is.na(out$hla_a) & !is.na(out$hla_b) & !is.na(out$hla_c))
  
  dat2 <- mapdf_dat %>% select(-hla_a, -hla_b, -hla_c) %>% 
    mutate(hla_a = "1/1", b1 = 1, b2 = 2, C1 = 1, C2 = 1)
  out2 <- map_score_df(dat2)
  expect_true(all(c("hla_a", "hla_b", "hla_c") %in% names(out2)))
  expect_true(!is.na(out2$hla_a) & !is.na(out2$hla_b) & !is.na(out2$hla_c))
  
  dat3 <- mapdf_dat %>% select(-hla_a, -hla_b, -hla_c) %>% 
    mutate(hla_a = NA, b1 = 1, b2 = "", C1 = NA, C2 = 1)
  out3 <- map_score_df(dat3)
  expect_true(all(c("hla_a", "hla_b", "hla_c") %in% names(out3)))
  expect_true(is.na(out3$hla_a) & is.na(out3$hla_b) & is.na(out3$hla_c))
})


context("Raefi and Krieger score functions")


test_that("Krieger score example works", {
          dat <- structure(list(kir_2DL1 = 1L, kir_2DL2 = 1L, kir_2DL3 = "NEG", 
                                kir_3DL1 = "NEG", kir_3DL2 = 1L, kir_2DS1 = 1L, kir_2DS2 = "NEG", 
                                kir_2DS4 = 1L, kir_2DS4N = "NEG", kir_2DS5 = 1L, kir_3DS1 = 1L, 
                                C_class = "C1/C1", B_class = "Bw6/Bw4-80T", A1 = "01:01", 
                                A2 = "01:01", B1 = "01:01", B2 = "01:01", C1 = "01:01", C2 = "01:01"), class = "data.frame", row.names = c(NA, 
                                                                                                                                           -1L))
          out <- score_krieger(dat)
          expect_equal(out$kirl_score, 1)
          expect_equal(out$imkir_score, 3)
          expect_equal(out$wkir_score, 2.78)
          })

test_that("Raefi score example works", {
  dat <- structure(list(kir_2DL1 = 1L, kir_2DL2 = 1L, kir_2DL3 = "NEG", 
                        kir_3DL1 = "NEG", kir_3DL2 = 1L, kir_2DS1 = 1L, kir_2DS2 = "NEG", 
                        kir_2DS4 = 1L, kir_2DS4N = "NEG", kir_2DS5 = 1L, kir_3DS1 = 1L, 
                        C_class = "C1/C1", B_class = "Bw6/Bw4-80T", A1 = "01:01", 
                        A2 = "01:01", B1 = "01:01", B2 = "01:01", C1 = "01:01", C2 = "01:01"), class = "data.frame", row.names = c(NA, 
                                                                                                                                   -1L))
  out <- score_rafei(dat)
  expect_equal(out$rafei_inh_kl_matches_2cat, "<3")
  expect_equal(out$rafei_act_kl_matches_2cat, ">=1")
  expect_equal(out$rafei_inact_kl_matches_2cat, "fav")
})

test_that("Raefi and Krieger scores match those used in 17-02 analysis", {
  # Henning wrote the original functions for calculating Raefi and Krieger scores for 17-02, here are exported 20 random samples with their 
  # scores calculated (without any IDs) - run them through the CTUtools functions and check the scores match
  dat3 <- read.csv("scores_check.csv", stringsAsFactors = FALSE)
  # sample 11 in the datafile gives different results for 3DL1 - this is becuase the B_class is "unknown", and the ifelse statements in the
  # original function counted this as present ligand, but it should be missing ligand as per the new function
  dat3 <- dat3[-11, ]
  
  s_k <- score_krieger(dat3)
  s_r <- score_rafei(dat3)
  
  
  expect_equal(dat3$kirl_score, s_k$kirl_score)
  expect_equal(dat3$imkir_score, s_k$imkir_score)
  expect_equal(dat3$wkir_score, s_k$wkir_score)
  
  expect_equal(dat3$rafei_inh_kl_matches_2cat, s_r$rafei_inh_kl_matches_2cat)
  expect_equal(dat3$rafei_act_kl_matches_2cat, s_r$rafei_act_kl_matches_2cat)
  expect_equal(dat3$rafei_inact_kl_matches_2cat, s_r$rafei_inact_kl_matches_2cat)
  
  
  
  })

# 
# dat <- read.csv("tmpdat.csv", stringsAsFactors = FALSE)
# score_krieger(dat)
# score_raefi(dat)
# dat3 <- read.csv("tests/testthat/scores_check.csv", stringsAsFactors = FALSE)
# score_krieger(dat3)
# score_raefi(dat3)
