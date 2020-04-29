test_that("cds2kaks(., ., model = 'Li') calculates correct value", {
  cds1 <- Biostrings::DNAStringSet("ATGCAACATTGC")
  cds2 <- Biostrings::DNAStringSet("ATGTATTGC")
  cds1.cds2.kaks1 <- cds2kaks(cds1, cds2, model = "Li")
  cds1.cds2.kaks2 <- unlist(seqinr::kaks(dnastring2aln(cds2codonaln(cds1, cds2))))
  expect_equal(cds1.cds2.kaks1, cds1.cds2.kaks2)
})
