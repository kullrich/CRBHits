data(ath)
data(aly)
data(ath_aly_crbh)

test_that("rbh2kaks() outputs correct values", {
    rbh.pairs <- ath_aly_crbh
    rbh.pairs$crbh.pairs <- rbh.pairs$crbh.pairs[1:2, ]
    ath_aly_crbh.kaks <- rbh2kaks(rbh.pairs, ath, aly)
    MSA2dist::dnastring2kaks(c(ath[4], aly[365]), isMSA=FALSE)
    MSA2dist::dnastring2kaks(c(ath[5], aly[54]), isMSA=FALSE)
    expect_true(ath_aly_crbh.kaks[1,8]==MSA2dist::dnastring2kaks(
        c(ath[4], aly[365]), isMSA=FALSE)[5] &&
        ath_aly_crbh.kaks[2,8]==MSA2dist::dnastring2kaks(
        c(ath[5], aly[54]), isMSA=FALSE)[5])
})
