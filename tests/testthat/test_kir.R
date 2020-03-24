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

context("KIR3DL1_3DS1 assignment")

test_that("KIR3DL1_3DS1_assignment examples work", {
  dat <- data.frame(kir_3DL1 = c("001/095+004", "POS", "005"), kir_3DS1 = c("NEG", "NEG", "013/107"), stringsAsFactors = FALSE)
  expect_equal(KIR3DL1_3DS1_assignment(dat$kir_3DL1, dat$kir_3DS1), c("KIR3DL1-H", "unknown", "KIR3DL1-L"))
  # Same with factor levels
  dat <- data.frame(kir_3DL1 = c("001/095+004", "POS", "005"), kir_3DS1 = c("NEG", "NEG", "013/107"), stringsAsFactors = TRUE)
  expect_equal(KIR3DL1_3DS1_assignment(dat$kir_3DL1, dat$kir_3DS1), c("KIR3DL1-H", "unknown", "KIR3DL1-L"))
  
  
})

context("KIR3DL1_HLA_B inhibiting function")

test_that("KIR3DL1_HLA_B_inhibiting examples work", {
  dat <- data.frame(assignment_KIR3DL1 = c("KIR3DL1-L", "KIR3DL1-N", "unknown", "unknown"), HLA_B_group = c("Bw4 - 80T", "Bw4 - 80I", NA, "Bw4 - 80I"), stringsAsFactors = FALSE)
  expect_equal(KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group), c("Strong inhibiting", "Weak inhibiting/noninhibiting", NA, "unknown"))
  expect_equal(KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group, levels = 3), c("Strong inhibiting", "noninhibiting", NA, "unknown"))
  expect_equal(KIR3DL1_HLA_B_inhibiting(dat$assignment_KIR3DL1, dat$HLA_B_group, levels = 4), c("Strong inhibiting", "Educated, Uninhibited", NA, "unknown"))
  # With factor dataframe
  dat_factor <- data.frame(assignment_KIR3DL1 = c("KIR3DL1-L", "KIR3DL1-N", "unknown", "unknown"), HLA_B_group = c("Bw4 - 80T", "Bw4 - 80I", NA, "Bw4 - 80I"), stringsAsFactors = TRUE)
  expect_equal(KIR3DL1_HLA_B_inhibiting(dat_factor$assignment_KIR3DL1, dat_factor$HLA_B_group), c("Strong inhibiting", "Weak inhibiting/noninhibiting", NA, "unknown"))
  
  })
