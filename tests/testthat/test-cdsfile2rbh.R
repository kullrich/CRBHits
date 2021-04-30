test_that("cdsfile2rbh() outputs correct number of conditional reciprocal best hit pairs", {
  #CRBHits_root <- system.file(package = "CRBHits")
  #LastTempDir <- tempdir()
  #system(paste0("unzip -o ", CRBHits_root, "/extdata/last-1219.zip -d ", LastTempDir))
  #system(paste0("cd ", LastTempDir, "/last-1219/; make"))
  tmp.paths <- make_vignette()
  athfile <- system.file("fasta/ath.cds.fasta.gz", package = "CRBHits")
  alyfile <- system.file("fasta/aly.cds.fasta.gz", package = "CRBHits")
  #ath_aly_crbh <- cdsfile2rbh(athfile, alyfile, lastpath = paste0(LastTempDir, "/last-1219/bin/"))

  ath_aly_crbh <- cdsfile2rbh(athfile, alyfile, lastpath = tmp.paths[1])
  expect_equal(dim(ath_aly_crbh$crbh.pairs)[1], 799)
})
