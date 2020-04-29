test_that("cds2codonaln() outputs correct alignment", {
  cds1 <- Biostrings::DNAString("ATGCAACATTGC")
  cds2 <- Biostrings::DNAString("ATGCATTGC")
  cds1.cds2.codonaln <- cds2codonaln(cds1, cds2)
  expect_true(
    as.character(cds1.cds2.codonaln[1])=="ATGCAACATTGC" &&
    as.character(cds1.cds2.codonaln[2])=="ATG---CATTGC"
  )
})
