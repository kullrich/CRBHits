data(ath)
data(aly)

test_that("cds2rbh() outputs correct number of conditional reciprocal best hit pairs", {
    #CRBHits_root <- system.file(package = "CRBHits")
    #LastTempDir <- tempdir()
    #system(paste0("unzip -o ", CRBHits_root, "/extdata/last-1256.zip -d ",
    #LastTempDir))
    #system(paste0("cd ", LastTempDir, "/last-1256/; make"))
    tmp.paths <- make_vignette()
    #ath_aly_crbh <- cds2rbh(ath, aly, lastpath=paste0(LastTempDir,
    #"/last-1256/bin/"))
    ath_aly_crbh <- cds2rbh(ath, aly, lastpath=tmp.paths[1])
    expect_equal(dim(ath_aly_crbh$crbh.pairs)[1], 211)
})
