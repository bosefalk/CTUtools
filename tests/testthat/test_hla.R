context("HLA_C classification")

test_that("HLA-C classification examples work", {
  expect_equal(HLA_C_classification("01:02", "01:AWFCH"), "C1/C1")
  dat <- data.frame(C1 = c("01:02", "02:03", "04:10"), C2 = c("01:AWFCH", "07:59", "05:50"), stringsAsFactors = FALSE)
  expect_equal(HLA_C_classification(dat$C1, dat$C2), c("C1/C1", "C1/C2", "C2/C2"))
  dat_odd <- data.frame(C1 = c("01:02", "01:02", "01:02"), C2 = c(NA, "", "07:02:01"), stringsAsFactors = FALSE)
  expect_equal(HLA_C_classification(dat_odd$C1, dat_odd$C2), c(NA_character_, NA_character_, NA_character_))
})

test_that("HLA-C classification works with factor columns", {
  dat <- data.frame(C1 = as.factor(c("01:02", "02:03", "04:10")), C2 = as.factor(c("01:AWFCH", "07:59", "05:50")))
  expect_equal(HLA_C_classification(dat$C1, dat$C2), c("C1/C1", "C1/C2", "C2/C2"))
  
})

context("HLA_B classification")

test_that("HLA-B classification examples work", {
  expect_equal(HLA_B_classification("07:02", "08:02"), "Bw4 - 80T")
  dat <- data.frame(B1 = c("07:02", "07:02", "07:02"), B2 = c("07:36", "39:15", "08:02"))
  expect_equal(HLA_B_classification(dat$B1, dat$B2), c("Bw4 - 80I", "Bw6", "Bw4 - 80T"))
  dat_odd <- data.frame(B1 = c("07:02", "07:02", "07:02"), B2 = c(NA, "", "08:47"))
  expect_equal(HLA_B_classification(dat_odd$B1, dat_odd$B2), c(NA_character_, NA_character_, NA_character_))
})

test_that("HLA-B classification works with factor columns", {
  dat <- data.frame(B1 = as.factor(c("07:02", "07:02", "07:02")), B2 = as.factor(c("07:36", "39:15", "08:02")))
  expect_equal(HLA_B_classification(dat$B1, dat$B2), c("Bw4 - 80I", "Bw6", "Bw4 - 80T"))
  
})

context("HLA generic classification")


test_that("HLA_Classification examples work", {
  dat <- data.frame(HLA_C = c("01:02", "01:AWFCH"), HLA_B = c("07:02", NA),stringsAsFactors = FALSE)
  HLA_C_class <- HLA_C_class_load()
  expect_equal(HLA_Classification(dat$HLA_C, HLA_C_class), c("C1", "C1"))
  HLA_B_class <- HLA_B_class_load()
  expect_equal(HLA_Classification(dat$HLA_B, HLA_B_class), c("Bw6", NA))
  
  dat_factor <- data.frame(HLA_C = c("01:02", "01:AWFCH"), HLA_B = c("07:02", NA),stringsAsFactors = TRUE)
  expect_equal(HLA_Classification(dat_factor$HLA_C, HLA_C_class), c("C1", "C1"))
  expect_equal(HLA_Classification(dat_factor$HLA_B, HLA_B_class), c("Bw6", NA))
  
})


