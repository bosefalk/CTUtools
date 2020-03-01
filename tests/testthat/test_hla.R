context("HLA_C classification")

test_that("HLA-C classification examples work", {
  expect_equal(HLA_C_classification("01:02", "01:AWFCH"), "C1/C1")
  dat <- data.frame(C1 = c("01:02", "02:03", "04:10"), C2 = c("01:AWFCH", "07:59", "05:50"), stringsAsFactors = FALSE)
  expect_equal(HLA_C_classification(dat$C1, dat$C2), c("C1/C1", "C1/C2", "C2/C2"))
  dat_odd <- data.frame(C1 = c("01:02", "01:02", "01:02"), C2 = c(NA, "", "07:02:01"), stringsAsFactors = FALSE)
  expect_equal(HLA_C_classification(dat_odd$C1, dat_odd$C2), c(NA_character_, NA_character_, NA_character_))
})

