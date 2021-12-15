data(ath)
data(aly)
data(ath_aly_crbh)

test_that("rbh2kaks() outputs correct values", {
    rbh.pairs <- ath_aly_crbh
    rbh.pairs$crbh.pairs <- rbh.pairs$crbh.pairs[1:2, ]
    ath_aly_crbh.kaks <- rbh2kaks(rbh.pairs, ath, aly)
    cds2kaks(ath[4], aly[365])
    cds2kaks(ath[5], aly[54])
    expect_true(ath_aly_crbh.kaks[1,4]==cds2kaks(ath[4],
        aly[365])[1] && ath_aly_crbh.kaks[2,4]==cds2kaks(ath[5], aly[54])[1])
})
