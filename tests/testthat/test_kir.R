context("KIR Gene Copy Number")

test_that("KIR_det_GCN examples work", {
  expect_equal(KIR_det_GCN("004+003|006+010"), "2")
  expect_equal(KIR_det_GCN("004+003|006+010", return_numeric = TRUE), 2)
  expect_equal(KIR_det_GCN("POS"), "x")
  expect_equal(KIR_det_GCN("POS", return_numeric = TRUE), NA_real_)
  
  dat <- data.frame(kircol = c("004+003", "001", "008+008+008", "POS", "NEG", "", NA, "001|001+002"))
  expect_equal(KIR_det_GCN(dat$kircol), 
               c("2", "1", "3", "x", "0", NA_character_, NA_character_, "1|2"))
  expect_equal(KIR_det_GCN(dat$kircol, return_numeric = TRUE), 
               c(2, 1, 3, NA_real_, 0, NA_real_, NA_real_, 1.5))
  
})

test_that("KIR_det_GCN works when passed factor column", {
  dat <- data.frame(kircol = as.factor(c("004+003", "001", "008+008+008", "POS", "NEG", "", NA, "001|001+002")))
  expect_equal(KIR_det_GCN(dat$kircol), 
               c("2", "1", "3", "x", "0", NA_character_, NA_character_, "1|2"))
  
})
