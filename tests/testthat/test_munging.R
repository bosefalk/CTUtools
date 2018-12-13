context("munging functions examples")

test_that("shortern_allele examples work", {
          expect_equal(shorten_allele("01:02:01:03"), "01:02")
          expect_equal(shorten_allele("03:ADJRE"), "03:ADJRE")
          expect_equal(shorten_allele("01:02:01:03", fields = 3), "01:02:01")

})


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


