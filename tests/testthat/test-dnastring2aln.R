test_that("dnastring2aln() outputs alignment", {
  cds1 <- Biostrings::DNAString("ATGCAACATTGC")
  cds2 <- Biostrings::DNAString("ATGCATTGC")
  cds1.cds2.codonaln <- cds2codonaln(cds1, cds2)
  dnaaln1 <- dnastring2aln(cds1.cds2.codonaln)
  dnaaln2 <- seqinr::as.alignment(nb = 2, nam = c("cds1", "cds2"),
                                                  seq = c("atgcaacattgc", "atg---cattgc"), NA)
  expect_equal(dnaaln1, dnaaln2)
})
