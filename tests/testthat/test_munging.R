context("shorten_allele munging")

test_that("shortern_allele examples work", {
  dat <- data.frame(HLA_A1 = c("01:02:01:03", "03:ADJRE", "01:02:01:01G"), stringsAsFactors = FALSE)
  expect_equal(shorten_allele(dat$HLA_A1), c("01:02", "03:ADJRE", "01:02"))
  expect_equal(shorten_allele(dat$HLA_A1, fields = 3), c("01:02:01", "03:ADJRE", "01:02:01"))
  dat2 <- data.frame(HLA_A1 = c("01:02:01:03", "03:ADJRE", NA, "some_random_string", 2, ""), stringsAsFactors = FALSE)
  expect_equal(shorten_allele(dat2$HLA_A1), c("01:02", "03:ADJRE", NA, "some_random_string", "2", NA))
})

test_that("shortern_allele examples work when allele column is factor", {
  dat <- data.frame(HLA_A1 = c("01:02:01:03", "03:ADJRE", "01:02:01:01G"), stringsAsFactors = TRUE)
  expect_equal(shorten_allele(dat$HLA_A1), c("01:02", "03:ADJRE", "01:02"))
  expect_equal(shorten_allele(dat$HLA_A1, fields = 3), c("01:02:01", "03:ADJRE", "01:02:01"))
  dat2 <- data.frame(HLA_A1 = c("01:02:01:03", "03:ADJRE", NA, "some_random_string", 2, ""), stringsAsFactors = TRUE)
  expect_equal(shorten_allele(dat2$HLA_A1), c("01:02", "03:ADJRE", NA, "some_random_string", "2", NA))
})

test_that("shortern_allele works when input only one string", {
  expect_equal(shorten_allele("01:02:01:03"), "01:02")
})

test_that("shortern_allele is ok with dplyr::rowwise() %>% mutate for legacy purposes", {
  dat <- data.frame(HLA_A1 = c("01:02:01:03", "03:ADJRE", "01:02:01:01G"), stringsAsFactors = FALSE)
  dat_mut <- dat %>% rowwise %>% mutate(HLA_A1 = shorten_allele(HLA_A1))
  expect_equal(dat_mut$HLA_A1, c("01:02", "03:ADJRE", "01:02"))
})

context("first_KIR_field munging")


test_that("first_KIR_field examples work", {
  expect_equal(KIR_first_field("0010101/0020102+0020103/0020104|0030105/0030106/0030107"), "001/002+002|003")
  expect_equal(KIR_first_field("001/002"), "001/002")
  expect_equal(KIR_first_field("0010203+0020304|NEG"), "001+002|NEG")
  expect_equal(KIR_first_field("POS"), "POS")

})


context("munging input verification errors")

test_that("shortern_allele input fields", {
  expect_error(shorten_allele("01:02:01:03", fields = 0), "shortern_allele: fields can only be set to 1,2,3 or 4")
  expect_error(shorten_allele("01:02:01:03", fields = 1), NA)
  expect_error(shorten_allele("01:02:01:03", fields = 2), NA)
  expect_error(shorten_allele("01:02:01:03", fields = 5), "shortern_allele: fields can only be set to 1,2,3 or 4")
  expect_error(shorten_allele("01:02:01:03", fields = NA), "shortern_allele: fields can only be set to 1,2,3 or 4")
  expect_error(shorten_allele("01:02:01:03", fields = "2"), NA)
  expect_error(shorten_allele("01:02:01:03", fields = "two"), "shortern_allele: fields can only be set to 1,2,3 or 4")
})

context("KIR first field")

test_that("KIR_first_field examples work", {
  dat <- data.frame(kirstring = c("0010101/0020102+0020103/0020104|0030105/0030106/0030107", "001/002", "NEG", "0010203+0020304|NEG", "POS", "", NA), stringsAsFactors = FALSE)
  expect_equal(KIR_first_field(dat$kirstring), c("001/002+002|003", "001/002", "NEG", "001+002|NEG", "POS", NA, NA))
  # Check with factor levels
  dat_factor <- data.frame(kirstring = c("0010101/0020102+0020103/0020104|0030105/0030106/0030107", "001/002", "NEG", "0010203+0020304|NEG", "POS", "", NA), stringsAsFactors = TRUE)
  expect_equal(KIR_first_field(dat_factor$kirstring), c("001/002+002|003", "001/002", "NEG", "001+002|NEG", "POS", NA, NA))
  
  
})

