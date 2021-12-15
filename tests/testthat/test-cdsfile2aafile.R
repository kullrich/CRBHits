data(ath)

test_that("cdsfile2aafile() outputs AAStringSet", {
    cdsfile <- system.file("fasta/ath.cds.fasta.gz", package="CRBHits")
    aafile <- tempfile()
    cdsfile2aafile(cdsfile, aafile)
    aa <- Biostrings::readAAStringSet(aafile)
    expect_true(class(aa)=="AAStringSet")
})
