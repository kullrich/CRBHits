data(ath)

test_that("cds2aa() outputs AAStringSet", {
  expect_true(class(cds2aa(ath)) == "AAStringSet")
})
