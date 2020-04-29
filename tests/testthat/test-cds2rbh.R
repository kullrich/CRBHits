data(ath)
data(aly)

test_that("cds2rbh() outputs correct number of conditional reciprocal best hit pairs", {
  CRBHits_root <- system.file(package = "CRBHits")
  LastTempDir <- tempdir()
  system(paste0("unzip -o ", CRBHits_root, "/extdata/last-1060.zip -d ", LastTempDir))
  system(paste0("cd ", LastTempDir, "/last-1060/; make"))
  ath_aly_crbh <- cds2rbh(ath, aly, lastpath = paste0(LastTempDir, "/last-1060/src/"))
  expect_equal(dim(ath_aly_crbh$crbh.pairs)[1], 798)
})
