data(ath)
data(aly)
data(ath_aly_crbh)

test_that("rbh2kaks() outputs correct values", {
  rbh.pairs <- ath_aly_crbh$crbh.pairs[1:2,]
  ath_aly_crbh.kaks <- rbh2kaks(rbh.pairs, ath, aly)
  cds2kaks(ath[1], aly[282])
  cds2kaks(ath[2], aly[280])
  expect_true(ath_aly_crbh.kaks[1,3] == cds2kaks(ath[1], aly[282])[1] && ath_aly_crbh.kaks[2,3] == cds2kaks(ath[2], aly[280])[1])
})
